from io import StringIO
from itertools import zip_longest
from copy import deepcopy
from collections import defaultdict

import pandas as pd

from Bio.Seq import Seq
from Bio.Blast import Applications, NCBIXML
from Bio import SeqIO
from Bio import pairwise2

from pyblast import utils

class BlocBlastException(Exception) :
    pass

class BCLine() :

    clines = {
        "blastp" : Applications.NcbiblastpCommandline,
        "blastn" : Applications.NcbiblastnCommandline,
        "blastx" : Applications.NcbiblastxCommandline,
        "tblastn" : Applications.NcbitblastnCommandline,
        "tblastx" : Applications.NcbitblastxCommandline,
        "psiblast" : Applications.NcbipsiblastCommandline,
        "pstblastn" : Applications.NcbirpstblastnCommandline,
        "rpstblastn" : Applications.NcbirpstblastnCommandline,
        "deltablast" : Applications.NcbideltablastCommandline
    }

    def __init__(self, bkind, ** kwargs) :

        self.query = kwargs.get("query", None)
        if not self.query : 
            raise BlocBlastException("Query not filled")

        if bkind not in BCLine.clines :
            raise BlocBlastException("Bkind not found in available blast executable : %s" %(" ".join(BCLine.clines)))

        self.cline = BCLine.clines[bkind](** kwargs)

    def __setitem__(self, name, value) :
        for parameter in self.cline.parameters :
            if name in parameter.names :
                parameter.value = value
                return

        raise BlocBlastException("Parameter not found : '%s'" %(name))

    def __str__(self) :
        return str(self.cline)

    @staticmethod
    def grouper_fdata(fdata, gsize) :
        # based on itertools recipy
        # https://docs.python.org/3.6/library/itertools.html#itertools-recipes

        args = [iter(fdata)] * gsize
        groups = zip_longest(*args, fillvalue=None)
        for group in groups :
            yield [seq for seq in group if seq is not None]

    @utils.ftiming
    def run_single_core(self, funcname="Blast command line (single)", postfun=None, quiet=False) :
        if not quiet : print (self)
        res = self.cline()
        return postfun(res) if postfun else res

    @staticmethod
    def run_single_core_tmp(tmpbcline, postfun=None, quiet=False) :
        # since we're not able to pickle a BioPython commande line
        # we're juste gonna execute it as a string
        
        if not quiet : print (str(tmpbcline))
        res = (utils.run_cmdline(str(tmpbcline), funcname="Blast command line (ncore)", quiet=quiet), None)
        return postfun(res) if postfun else res

    def run_chunked(self, ncore, chunksize=200, postfun=None, quiet=False) :
        fdata = SeqIO.parse(self.query, "fasta")
        fdata = BCLine.grouper_fdata(fdata, chunksize)
        queries = (utils.TMPFasta(group, quiet=quiet) for group in fdata)

        fun = BCLine.run_single_core_tmp
        args = [utils.FunArgs(fun, TMPBCLine(self, query), postfun=postfun, quiet=quiet) for query in queries]
        
        return utils.mproc(args, ncore)

    def run(self, ncore=1, chunksize=200, postfun=None, quiet=False) :
        if ncore == 1 :
            return self.run_single_core(postfun=postfun, quiet=quiet)
        else :
            return self.run_chunked(ncore, chunksize, postfun=postfun, quiet=quiet)

    @staticmethod
    def mbdb(db, ** kwargs) :
        # Applications.NcbimakeblastdbCommandline not found
        kwargs["in"] = db
        cmd = "makeblastdb %s" %(" ".join("-%s %s" %(key, value) for key, value in kwargs.items()))
        utils.run_cmdline(cmd)

    @staticmethod
    def clean_fasta(fname, astmp=True) :
        # Change header to match correct id with _
        # Replace ? to N
        # Return a new fasta (tmp)

        def clean_seqrecord(seqrecord) :
            name = seqrecord.description.replace(" ", "_")
            seqrecord.id = seqrecord.name = seqrecord.description = name
            sequence = str(seqrecord.seq).upper().replace("?", "N")
            seqrecord.seq = Seq(sequence)
            return seqrecord

        fdata = SeqIO.parse(fname, "fasta")
        fdata = (clean_seqrecord(sr) for sr in fdata)
        return utils.TMPFasta(fdata) if astmp else fdata

class TMPBCLine() :

    # A single class used when we multiprocess blast
    # We cannot pickle BioPython command line
    # And we have to keep the TMPFasta instance alive otherwise the file is removed

    def __init__(self, bcline, tmpfasta) :
        bcline = deepcopy(bcline)
        bcline["query"] = str(tmpfasta)

        self.tmpfasta = tmpfasta
        self.cmdline = str(bcline)

    def __str__(self) :
        return self.cmdline

class BCLine6(BCLine) :

    def __init__(self, bkind, ** kwargs) :

        kwargs, self.spec = self.add_specifiers(kwargs)
        super(BCLine6, self).__init__(bkind, ** kwargs)

    def add_specifiers(self, kwargs) :
        current = kwargs.get("outfmt", "").split(" ")
        current = [value for value in current if value]

        toadd = ["qseqid", "sseqid", "pident", "nident", "length", "slen", 
        "qlen", "qstart", "qend", "sstart", "send", "positive"]
        used = list(set(current) | set(toadd))

        kwargs["outfmt"] = "'%s'" %("6 " + " ".join(used))
        return kwargs, used

    def treat_res(self, res) :

        if isinstance(res, list) :
            # multiprocessing res
            return pd.concat([self.treat_res(r) for r in res])

        out, err = res
        if err : print (err)

        out = StringIO(out)
        df = pd.read_csv(out, sep="\t", names=self.spec)

        # check if a column is full of na
        # may happend if one gave a wrong outfmt
        for column in df.columns :
            if df[column].isnull().all() and not df.empty :
                print ("""WARNING : An output column is full of nan :
                    This could be a result of a wrong outfmt key,
                    and may break result outputs\n""")

        df["qlength"] = abs(df["qend"] - df["qstart"]) + 1
        df["slength"] = abs(df["send"] - df["sstart"]) + 1

        prc = [("qide", "nident", "qlen"), ("qcov", "qlength", "qlen"), ("qpos", "positive", "qlen"), 
               ("side", "nident", "slen"), ("scov", "slength", "slen"), ("spos", "positive", "slen")]

        for name, divident, divisor in prc :
            df[name] = df[divident] * 100 / df[divisor]

        return df

    @staticmethod
    def get_best(df, query="qseqid", column="qpos", nmatches=10, use_max=True) :
        if use_max :
            return BCLine6.get_best_max(df, query, column, nmatches)
        else :
            return BCLine6.get_best_min(df, query, column, nmatches)

    @staticmethod
    def get_best_max(df, query="qseqid", column="qpos", nmatches=10) :
        fun = (lambda serie : serie.nlargest(nmatches).min()) if nmatches > 1 else max
        mask = df.groupby(query)[column].transform(fun)
        return df[df[column] >= mask]

    @staticmethod
    def get_best_min(df, query="qseqid", column="qpos", nmatches=10) :
        fun = (lambda serie : serie.smallest(nmatches).max()) if nmatches > 1 else min
        mask = df.groupby(query)[column].transform(fun)
        return df[df[column] <= mask]    

    @staticmethod
    def query_mdbs(bkind, query, dbs, best=False, chunksize=200, ncore=1,
        mcolumn="qpos", nmatches=1, ** kwargs) :

        res = []
        for db in dbs :
            df = BCLine6(bkind, query=query, subject=db, ** kwargs).run(ncore=ncore, chunksize=chunksize)
            if df.empty : continue
            df["db"] = db
            res.append(df)

        if not res : return pd.DataFrame()
        res = pd.concat(res)

        return BCLine6.get_best(res, column=mcolumn, nmatches=nmatches) if best else res

    def run(self, * args,  ** kwargs) :
        res = super(BCLine6, self).run(* args, ** kwargs)
        return self.treat_res(res)

class BCLineHSPFuse(BCLine) :

    def __init__(self, bkind, ** kwargs) :

        kwargs = self.add_specifiers(kwargs)
        super(BCLineHSPFuse, self).__init__(bkind, ** kwargs)

    def add_specifiers(self, kwargs) :
        current = kwargs.get("outfmt", "").split(" ")
        current = current[1:] if len(current) > 1 else []

        kwargs["outfmt"] = "'%s'" %("5 " + " ".join(current))
        return kwargs

    @staticmethod
    def parse_xml(data) :
        f = StringIO(data)
        try : return [record for record in NCBIXML.parse(f)]
        except ValueError : return [] # empty results

    @staticmethod
    def iter_hsp(records) :
        for record in records :
            qname, qlength = record.query, record.query_length
            for aln in record.alignments :
                sname, slength = aln.hit_def, aln.length
                for hsp in aln.hsps :
                    yield (qname, qlength), (sname, slength), hsp

    @staticmethod
    def get_hsp_ident(hsp, query=True) :
        start = hsp.query_start if query else hsp.sbjct_start
        end = hsp.query_end if query else hsp.sbjct_end
        sequence = hsp.query if query else hsp.sbjct
        reverse = start > end

        positives, ident = set(), set()
        psequence = -1

        for idx, letter in enumerate(sequence) :
            letter_match = hsp.match[idx]

            if letter_match != ' ' :
                position = start - psequence if reverse else psequence + start
                positives.add(position)
                if letter_match != "+" : ident.add(position)

            psequence = psequence + 1 if letter != "-" else psequence

        return positives, ident

    @staticmethod
    def get_sum_sim_all(records) :
        qlength, slength = {}, {}
        data = defaultdict(lambda : defaultdict(list))

        for query, subject, aln in BCLineHSPFuse.iter_hsp(records) :
            qname, sname = query[0], subject[0]
            qlen, slen = query[1], subject[1]
            qlength[qname], qlength[sname] = qlen, slen

            qposit, qident = BCLineHSPFuse.get_hsp_ident(aln, query=True)
            sposit, sident = BCLineHSPFuse.get_hsp_ident(aln, query=False)
            positions = {"qide" : qident, "qpos" : qposit, "side" : sident, "spos" : sposit}
            data[qname][sname].append(positions)

        df = []
        kinds = ["qide", "qpos", "side", "spos"]

        for qname, subject in data.items() :
            for sname, all_positions in subject.items() :

                hsp_count = len(all_positions)

                if len(all_positions) > 1 :
                    intersection = {kind + "_inter" : len(set.intersection(* [positions[kind]
                    for positions in all_positions])) for kind in kinds}

                else :
                    intersection = {kind + "_inter" : 0 for kind in kinds}

                union = {kind : len(set.union(* [positions[kind]
                for positions in all_positions])) for kind in kinds}

                info = {"qseqid" : qname, "sseqid" : sname, "#hsp" : hsp_count, ** union, ** intersection}
                df.append(info)

        df = pd.DataFrame(df)
        if df.empty : return df

        df["qlen"], df["slen"] = df["qseqid"].map(qlength), df["sseqid"].map(qlength)
        df["qide_prc"] = df["qide"] * 100 / df["qlen"]
        df["qpos_prc"] = df["qpos"] * 100 / df["qlen"]
        df["side_prc"] = df["side"] * 100 / df["slen"]
        df["spos_prc"] = df["spos"] * 100 / df["slen"]

        df = df[sorted(df.columns)]
        return df

    @staticmethod
    def post_treat(res) :
        out, err = res
        records = BCLineHSPFuse.parse_xml(out)
        return BCLineHSPFuse.get_sum_sim_all(records)

    def treat_res(self, res) :

        if isinstance(res, list) :
            # multiprocessing res
            return pd.concat(res)

        return res

    def run(self, * args,  ** kwargs) :
        res = super(BCLineHSPFuse, self).run(* args, postfun=BCLineHSPFuse.post_treat, ** kwargs)
        return self.treat_res(res)

class GlobalAlignment() :

    def __init__(self, inputs) :
        # inputs are an iterable of pairs of SeqRecords
        self.inputs = inputs

    @staticmethod
    def get_alignfun(name="globalxx") :
        return getattr(pairwise2.align, name)

    @staticmethod
    def aln_identity(query, subject, gapless=True) :
        if gapless :
            return sum(1 if b1 == b2 and b1 != "-" and b2 != "-" else 0 for b1, b2 in zip(query, subject))
        else :
            return sum(1 if b1 == b2 else 0 for b1, b2 in zip(query, subject))

    @staticmethod
    def run_pair(query, subject, * args, gfun="globalxx", ** kwargs) :
        gfun = GlobalAlignment.get_alignfun(gfun)
        alns = gfun(query.seq, subject.seq, * args, ** kwargs)
        return query, subject, alns

    @staticmethod
    def run_pair_identity(query, subject, * args, gapless=True, addnames=True, ** kwargs) :
        query, subject, alns = GlobalAlignment.run_pair(query, subject, * args, ** kwargs)
        values = []

        for aln in alns :
            alnID = GlobalAlignment.aln_identity(aln[0], aln[1], gapless=gapless)
            values.append({"gID" : alnID, "gPID" : alnID / len(aln[0]), "gQPID" : alnID * 100 / len(query), 
            "gSPID" : alnID * 100 / len(subject), "gSCORE" : aln[2], "gLEN" : len(aln[0]),
            "qseqid" : query.id, "sseqid" : subject.id})

        return query, subject, pd.DataFrame(values)

    @staticmethod
    def run_pair_score(query, subject, * args, ** kwargs) :
        kwargs["one_alignment_only"] = True
        kwargs["score_only"] = False
        return GlobalAlignment.run_pair(query, subject, * args, ** kwargs)

    def run(self, * args, ncore=1, callback=None, ** kwargs) :
        callback = callback or GlobalAlignment.run_pair
        args = [utils.FunArgs(callback, query, subject, * args, ** kwargs) for query, subject in self.inputs]
        return utils.mproc(args, ncore=ncore)

class Local2Global() :

    def __init__(self, query, subject) :
        # Perform local and global alignment based on the best results
        self.query = query
        self.subject = subject

    def retrieve_seq(self, query=True, ids=None) :
        fname = self.query if query else self.subject
        seqs = SeqIO.parse(fname, "fasta")
        if ids : return (seq for seq in seqs if seq.id in ids)
        return seqs

    def run(self, bkind, gargs=(), gkwargs={}, bkwargs={}, match_fun=BCLine6.get_best, chunksize=200, ncore=1, quiet=False) :
        # Blast local alignments
        bc = BCLine6(bkind, query=self.query, db=self.subject, ** bkwargs)
        bres = bc.run(chunksize=chunksize, ncore=ncore, quiet=quiet)
        if match_fun : bres = match_fun(bres)

        # retrieve subject and query sequences
        seq_que = {seq.id : seq for seq in self.retrieve_seq(query=True, ids=set(bres["qseqid"].unique()))}
        seq_sub = {seq.id : seq for seq in self.retrieve_seq(query=False, ids=set(bres["sseqid"].unique()))}

        # we run global alignments
        pairs = bres[["qseqid", "sseqid"]].drop_duplicates(["qseqid", "sseqid"])        
        seqs = [(seq_que[query], seq_sub[subject]) for query, subject in zip(pairs["qseqid"], pairs["sseqid"])]
        
        print ("%i global alignments to process" %(len(seqs)))
        call = GlobalAlignment.run_pair_identity
        gres = GlobalAlignment(seqs).run(* gargs, ncore=ncore, callback=call, ** gkwargs)
        gres = pd.concat([res[2] for res in gres])

        return bres.merge(gres, on=["qseqid", "sseqid"], how="left")