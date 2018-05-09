###################################################################################################
# Author: Skyler Kuhn, NIH/NCI [C], Leidos Biomedical Research Inc.
# Web-scraper to ultimately convert GEO accession numbers to SRA Run IDs-- needed for sratoolkit
## Usage: python GEOtoSRArun.py
## Dependencies: bs4, requests (if needed use: pip install)
###################################################################################################

from __future__ import print_function
from bs4 import BeautifulSoup as bs
import requests


def convertGEO2SRA(geo):
    """
    Writes GEO Accession Number, SRArunID, desc to an output file
    :param geo:
    :yield SRArunID:
    """
    # Using the GEO ID to goto the correct Webpage
    request = requests.get('https://www.ncbi.nlm.nih.gov/sra/?term={}'.format(geo))
    content = request.content
    soup = bs(content, 'html.parser')

    # Parsing the html to find the SRA_run_id (used by sratoolkit)
    table = soup.findChildren('table')[0]
    rows = table.findChildren('tr')

    for row in rows:
        try:
            SRA_run_id = row.findChildren('td')[0]
            yield SRA_run_id.string
        except IndexError:
            pass


def geo_accensions(filename):
    """
    Pulls the GeoID and Description from the given TSV file
    :param filename:
    :yield geo accession number, description:
    """
    fh = open(filename)
    for line in fh:
        linelist = line.strip().split("\t")
        geoid, desc = linelist[0], linelist[1]
        yield geoid, desc
    fh.close()


if __name__ == "__main__":
    outfh = open("GEO_SRA_DESC.tsv", "w")
    outSwarm = open("downloadfastq.swarm","w")
    outSwarm.write("#Swarm Usage: swarm -f downloadfastq.swarm --module sratoolkit --gres=lscratch:100 -g 32\n")
    filecount = 0
    for geoacc, description in geo_accensions("ChIPSeqBenchmarkingDatasets.txt"):
        for SRArun in convertGEO2SRA(geo=geoacc):
            filecount += 1
            print("Please wait... webpage {} is being scraped".format(filecount))
            outfh.write("{}\t{}\t{}\n".format(geoacc, SRArun, description))
            outSwarm.write("cd /scratch/ChIPSeqBenchmarking; fastq-dump --split-files --gzip {}\n".format(SRArun)) # change this later so program accepts target dump location

    outfh.close()
    outSwarm.close()
    print("Done!")
