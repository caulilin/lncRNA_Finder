lncRNA_Finder
=============

The open source code of pipeline for the identification of long non-coding RNAs

Source code of the pipeline-LncRNA_Finder

1, Copyright (c) 2013

University of Minnesota, All Rights Reserved. Authors: Lin Li, Nathan M. Springer, Gary J. Muehlbauer

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

   o Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

   o Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

   o Neither the name of the University of Minnesota nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LIN LI, NATHAN M. SPRINGER, or GARY J. MUEHLBAUER (OR UNIVERSITY OF MINNESOTA) BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR  TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

DESCRIPTION: LncRNA_Finder enables the discovery of long noncoding RNAs using native sequence fasta files. This script essentially uses the results from external alignment programs and performs long noncoding RNA filtering via a set of specified parameters.

CITATION: LncRNA_Finder can be cited as: Li L, Eichten SR, Shimizu R, Petsch K, Yeh CT, Wu W, Chettoor AM,  Givan SA, Cole RA, Fowler JE, Evans MMS, Scanlon MJ, Yu JM, Schnable PS, Timmermans MCP, Springer NM, and Muehlbauer GJ: Genome-wide discovery and characterization of maize long non-coding RNAs. Genome Biology 2014, 15:R40.

2, Prerequisite of external softwares

Need to install ncbi_blast standalone package, bowtie and cpc in your local environment correctly and set the running path of external programs at the begining of the pipeline

3, Usage

Perl  LncRNA_Finder.pl -i \<transcript.fasta\> -p \<protein.fasta\> -k \<housekeeping.fasta\> -s \<smallRNA.fasta\> -o \<output prefix\> [-t <# of thread>] [-r <minimum lncRNA length>] [-f <maximum ORF length>] [-m <# of mismatch>] [-e <E-value of alignment>]

Options:

   -i \<transcript.fasta\>

   -p \<protein.fasta\>

   -k \<housekeeping.fasta\>

   -s \<smallRNA.fasta\>

   -o \<output prefix\>

   -h help

   -t <int> number of thread for the computation || default=4

   -r <int> minimum lncRNA length || default=200

   -f <int> maximum potential ORF length of lncRNAs || default=100

   -m <int> number of mismatch in the alignment with smallRNA || default=0

   -e E-value of the alignment against protein database || default=1.0e-3

