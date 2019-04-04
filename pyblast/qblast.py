# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2019-04-03 16:04:11
# @Last modified by:   jsgounot
# @Last Modified time: 2019-04-04 10:55:04

import os, glob
import time, datetime
from threading import Thread

import pandas as pd

from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

from pyblast import BCLineHSPFuse
from pyblast import utils

"""

Online blast query

"""

class QBlast(object):

    def __init__(self, fasta, * args, xml_file=None, ** kwargs) :

        self.fasta = fasta

        if xml_file is not None :
            self.xml_file = xml_file

        else :
            with open(self.fasta) as f :
                kwargs["sequence"] = f.read()

            self.xml_file = self.query_ncbi(* args, ** kwargs)

    def query_ncbi(self, * args, ** kwargs) :

        res = NCBIWWW.qblast(* args, ** kwargs)
        xml_file = utils.TMPFname()

        with open(str(xml_file), "w") as f :
            f.write(res.read())

        res.close()
        return xml_file

    @staticmethod
    def parse_xml_file(filename) :
        with open(filename) as f :
            try : return [record for record in NCBIXML.parse(f)]
            except ValueError : return [] # empty results

    def as_table(self) :
        records = QBlast.parse_xml_file(str(self.xml_file))
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
        records = QBlast.parse_xml_file(str(self.xml_file))
        return BCLineHSPFuse.get_sum_sim_all(records)

    @staticmethod
    def thread_qblast(fasta, results, * args, ** kwargs) :
        print ("%s : Run qblast for : %s" %(datetime.datetime.now(), fasta))
        results.append(QBlast(fasta, * args, ** kwargs))
        print ("%s : Finish blast for : %s" %(datetime.datetime.now(), fasta))

    @staticmethod
    def threads_qblast(fasta_list, * args, time_wait=10, ** kwargs) :
        # Run parallel qblast queries over ncbi using threads
        # https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=DeveloperInfo
        # time wait corresponds to the time beetween two parallele requests, if you don't want that
        # ncbi kicks your ass, do not lower this value

        results = []
        threads = []
        
        fasta_list = glob.glob(fasta_list) if isinstance(fasta_list, str) else fasta_list
        fun = lambda fasta, results : QBlast.thread_qblast(fasta, results, * args, ** kwargs)      

        for fasta in fasta_list :
            thread = Thread(target=fun, args=(fasta, results))
            thread.start()
            threads.append(thread)
            time.sleep(time_wait)

        for thread in threads :
            thread.join()

        return {res.fasta : res for res in results}