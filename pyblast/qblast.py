# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2019-04-03 16:04:11
# @Last modified by:   jsgounot
# @Last Modified time: 2019-04-08 09:56:37

import os, glob, shutil
import time, datetime

from itertools import chain
from threading import Thread

import pandas as pd

from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

from pyblast import BCLine, BCLineHSPFuse
from pyblast import utils

"""

Online blast query

"""

current_time = lambda : datetime.datetime.now().strftime("%H:%M:%S")

class Query(object) :

    def __init__(self, qbr_instance, job_name, program, database, fasta_str, ** kwargs) :
        self.qbr_instance = qbr_instance
        self.job_name = job_name

        self.program = program
        self.database = database
        self.fasta_str = fasta_str
        self.kwargs = kwargs

        if kwargs.get("format_type", "XML") != "XML" :
            raise ValueError("Only XML output are managed")

    def trigger(self) :
        print ("%s : Run qblast (%s - %s) - job name : %s" %(current_time(), self.program, self.database, self.job_name))
        res = NCBIWWW.qblast(self.program, self.database, self.fasta_str, ** self.kwargs)
        print ("%s : End qblast - job name : %s" %(current_time(), self.job_name))

        xml_file = utils.TMPFname()
        with open(str(xml_file), "w") as f :
            f.write(res.read())

        res.close()
        self.qbr_instance.append(xml_file)

class QueriesManager(object) :

    """

    Run multiple queries within threads

    """

    def __init__(self, queries=[]) :
        self.queries = queries

    def __add__(self, other) :
        if not isinstance(other, QueriesManager) :
            raise ValueError("Queries Manager can only be associated with other QM")

        new_queries = self.queries + other.queries
        return QueriesManager(new_queries)

    def run(self, time_wait=15) :
        # Run parallel qblast queries over ncbi using threads
        # https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=DeveloperInfo
        # time wait corresponds to the time beetween two parallele requests, if you don't want that
        # ncbi kicks your ass, do not lower this value

        threads = []
        fun = lambda query : query.trigger()

        for query in self.queries :
            thread = Thread(target=fun, args=(query,))
            thread.start()
            threads.append(thread)
            time.sleep(time_wait)

        for thread in threads :
            thread.join()

    # --------------------------------------------------
    # Create the Queries Manager from data (fasta)
    # --------------------------------------------------

    @staticmethod
    def records_to_str(records) :
        strs = []
        
        for record in records :
            identifier, seq = record.id, str(record.seq)
            strs.append(">%s" %(identifier))

            for i in range(0, len(seq), 60) :
                strs.append(seq[i:i+60])

        return "\n".join(strs)

    @staticmethod
    def extract_fasta_str(fasta, chunk_size=10) :
        # return the fasta sequence as a
        # string object

        fdata = list(SeqIO.parse(fasta, "fasta"))

        if chunk_size is None :
            return [QueriesManager.records_to_str(fdata)]

        return [QueriesManager.records_to_str(records)
            for records in BCLine.grouper_fdata(fdata, chunk_size)]

    @staticmethod
    def from_fasta(program, database, fasta, chunk_size=None, ** kwargs) :
        # Return the QBlastResult and the Queries manager based on qblast parameter and fasta chunk_size

        qbr = QBlastResult(fasta=fasta)
        queries = []
        subsequences = list(QueriesManager.extract_fasta_str(fasta, chunk_size))

        for idx, subsequence in enumerate(subsequences, start=1) :
            job_name = "%s (n = %i)" %(fasta, idx) if len(subsequences) > 1 else fasta
            queries.append(Query(qbr, job_name, program, database, subsequence, ** kwargs))

        return qbr, QueriesManager(queries)

    @staticmethod
    def from_multiple_fastas(program, database, fasta_list, chunk_size=None, ** kwargs) :
        qbrs, qms = zip(* [QueriesManager.from_fasta(program, database, fasta, chunk_size, ** kwargs)
            for fasta in fasta_list])

        qm = sum(qms, QueriesManager())
        qbrs = qbrs if len(qbrs) > 1 else next(iter(qbrs))

        return qbrs, qm

class QBlastResult(object):

    """

    Wrapper of NCBIWWW Qblast function 
    http://biopython.org/DIST/docs/api/Bio.Blast.NCBIWWW-module.html

    """

    def __init__(self, xml_files=None, fasta=None) :

        self.fasta = fasta

        if xml_files is not None :
            xml_files = glob.glob(xml_files) if isinstance(xml_files, str) else xml_files
            self.xml_files = xml_files

        else :
            self.xml_files = []

    def append(self, xml_file) :
        self.xml_files.append(xml_file)

    @staticmethod
    def fuse_xml(xml_files) :
        raise Exception()

    def save(self, outfile) :
        if len(self.xml_files) == 0 : raise ValueError("Empty QBlastResult instance")
        xml_file = self.xml_files[0] if len(self.xml_files) == 1 else fuse_xml(self.xml_files)
        shutil.copy_file(xml_file, outfile)

    @staticmethod
    def parse_xml_file(xml_file) :
        with open(str(xml_file)) as f :
            try : return [record for record in NCBIXML.parse(f)]
            except ValueError : return [] # empty results

    def parse_xml_files(self) :

        if not self.xml_files :
            raise ValueError("No xml files produced or provided")

        return chain(* (QBlastResult.parse_xml_file(xml_file)
            for xml_file in self.xml_files))

    def as_table(self) :
        records = self.parse_xml_files()
        rows = []

        for query, subject, hsp in BCLineHSPFuse.iter_hsp(records) :
            row = {}

            row["qname"], row["sname"] = query[0], subject[0]
            row["qlen"], row["slen"] = query[1], subject[1]
            row["length"], row["gap"] = hsp.align_length, hsp.gaps
            row["hsp_ident"], row["hsp_pos"] = hsp.identities, hsp.positives
            row["expect"] = hsp.expect

            rows.append(row)

        df = pd.DataFrame(rows)
        if df.empty : return df

        df["ppos"] = df["hsp_pos"] * 100 / df["length"]
        df["pident"] = df["hsp_ident"] * 100 / df["length"]

        df["qide_prc"] = df["hsp_ident"] * 100 / df["qlen"]
        df["qpos_prc"] = df["hsp_pos"] * 100 / df["qlen"]
        df["side_prc"] = df["hsp_ident"] * 100 / df["slen"]
        df["spos_prc"] = df["hsp_pos"] * 100 / df["slen"]

        return df

    def hsp_fuse(self) :
        records = self.parse_xml_files()
        return BCLineHSPFuse.get_sum_sim_all(records)


# --------------------------------------------------
# Called function
# --------------------------------------------------

def launch_qblasts(program, database, fasta_list, chunk_size=None, time_wait=15, ** kwargs) :
    # Run parallel or single qblast queries over ncbi using threads
    # https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=DeveloperInfo
    # time wait corresponds to the time beetween two parallele requests, if you don't want that
    # ncbi kicks your ass, do not lower this value

    fasta_list = glob.glob(fasta_list) if isinstance(fasta_list, str) else fasta_list
    qbr, qm = QueriesManager.from_multiple_fastas(program, database, fasta_list, chunk_size=chunk_size, ** kwargs)
    qm.run(time_wait)
    return qbr
