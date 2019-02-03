from io import StringIO
from itertools import zip_longest
from copy import deepcopy

import pandas as pd

from Bio.Seq import Seq
from Bio.Blast import Applications
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

    def grouper_fdata(self, fdata, gsize) :
        # based on itertools recipy
        # https://docs.python.org/3.6/library/itertools.html#itertools-recipes

        args = [iter(fdata)] * gsize
        groups = zip_longest(*args, fillvalue=None)
        for group in groups :
            yield [seq for seq in group if seq is not None]

    @utils.ftiming
    def run_single_core(self, funcname="Blast command line (single)") :
        print (self)
        return self.cline()

    @staticmethod
    def run_single_core_tmp(tmpbcline) :
        # since we're not able to pickle a BioPython commande line
        # we're juste gonna execute it as a string
        print (str(tmpbcline))
        return (utils.run_cmdline(str(tmpbcline), funcname="Blast command line (ncore)"), None)

    def run_chunked(self, ncore, chunksize=200) :
        fdata = SeqIO.parse(self.query, "fasta")
        fdata = self.grouper_fdata(fdata, chunksize)
        queries = (utils.TMPFasta(group) for group in fdata)

        fun = BCLine.run_single_core_tmp
        args = [utils.FunArgs(fun, TMPBCLine(self, query)) for query in queries]
        
        return utils.mproc(args, ncore)

    def run(self, ncore=1, chunksize=200) :
        if ncore == 1 :
            return self.run_single_core()
        else :
            return self.run_chunked(ncore, chunksize)

    @staticmethod
    def mbdb(db, ** kwargs) :
        # Applications.NcbimakeblastdbCommandline not found
        kwargs["in"] = db
        cmd = "makeblastdb %s" %(" ".join("-%s %s" %(key, value) for key, value in kwargs.items()))
        utils.run_cmdline(cmd)

    @staticmethod
    def clean_fasta_nucl(fname, minsize=None) :
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
        return utils.TMPFasta(fdata)


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
        current = current[1:] if len(current) > 1 else []

        toadd = ["qseqid", "sseqid", "pident", "nident", "length", "slen", "qlen", "qstart", "qend", "sstart", "send", "positive"]
        used = list(set(current) | set(toadd))

        kwargs["outfmt"] = "'%s'" %("6 " + " ".join(used))
        return kwargs, used

    def treat_res(self, res) :

        if isinstance(res, list) :
            # multiprocessing res
            return pd.concat([self.treat_res(r) for r in res])

        cols = list(self.spec)

        out, err = res
        if err : print (err)

        out = StringIO(out)
        df = pd.read_csv(out, sep="\t", names=cols)

        df["qlength"] = abs(df["qend"] - df["qstart"]) + 1
        df["slength"] = abs(df["send"] - df["sstart"]) + 1

        prc = [("qide", "nident", "qlen"), ("qcov", "qlength", "qlen"), ("qpos", "positive", "qlen"), 
               ("side", "nident", "slen"), ("scov", "slength", "slen"), ("spos", "positive", "slen")]

        for name, divident, divisor in prc :
            df[name] = df[divident] * 100 / df[divisor]

        return df

    @staticmethod
    def get_best(df, column="qpos", nmatches=10) :
        fun = (lambda serie : serie.nlargest(nmatches).min()) if nmatches > 1 else max
        mask = df.groupby("qseqid")[column].transform(fun)
        return df[df[column] >= mask] 

    def run(self, * args,  ** kwargs) :
        res = super(BCLine6, self).run(* args, ** kwargs)
        return self.treat_res(res)

    @staticmethod
    def query_mdb(bkind, query, dbs, best=False, chunksize=200, ncore=1,
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

    def run(self, bkind, gargs=(), gkwargs={}, bkwargs={}, match_fun=BCLine6.get_best, chunksize=200, ncore=1) :
        # Blast local alignments
        bc = BCLine6(bkind, query=self.query, db=self.subject, ** bkwargs)
        bres = bc.run(chunksize=chunksize, ncore=ncore)
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