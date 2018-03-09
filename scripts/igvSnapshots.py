##########################################################################################################
# ChIP-Seq Pipeline: MAC Peakcaller -> IGV snapshot
# Author: Skyler Kuhn (NIH/NCI) [C]
# CCR Collaborative Bioinformatics Resource
# Version 1.0.2, Excepted Release: 01/30/2017
# See readme.txt for more information
# USAGE:
#   python igvSnapshots.py
#    --n=5
#    --narrowPeak_file=CHIP_Thpok_Biotin_vs_Input_Thpok_peaks.narrowPeak
#    --treatmentBW_file=CHIP_Thpok_Biotin.R1.trim.not_blacklist_plus.sorted.mapq_gt_3.normalized.bw
#    --inputBW_file=Input_Thpok.R1.trim.not_blacklist_plus.sorted.mapq_gt_3.normalized.bw
#    --output_folder=SNAPTEST
##########################################################################################################

# Imports
from __future__ import print_function, division
import subprocess as sp
import sys
import time
import os


def check_args(all_args):
    """
    :param all_args: # (this is a list of all provided command-line arguements)
    :return: arg_dict # if 6 or 11arguments are not provided an Exception is raised
    TLDR: This function checks the provided command-line arguments and them uses regex to check to see if they are valid,
    if they are not an Exception is raised!
    """
    def parse_args(args):  # maybe look into using a decorator here, the more python-ic way
        """
        Generator that returns,
        :param args: (the list of args provided by sys.argv)
        :return:
        """
        for i in range(1, len(args)-1, 2):  # the zero-th index is the program name
            j = i + 1
            yield args[i], args[j]

    def invalid_usage_message():
        """
        :return: docstring error message
        TDLR: An error Message is returned if the wrong number of command-line args are provided
        """
        return """Failed to Provide Required Input Arguments:
            --n
            --narrowPeak_file
            --treatmentBW_file
            --inputBW_file
            --output_folder\n* Invalid Input Arguments provided *
            \nUsage:\npython igvSnapshots.py --n=5 --narrowPeak_file=narrow.narrowPeak --treatmentBW_file=treatment.bw --inputBW_file=inputfile.bw --output_folder=FolderName
            """

    if len(all_args) != 7 and len(all_args) != 13:  # maybe use an assert statement here
        raise Exception(invalid_usage_message())

    arg_dict = {}  # k:switch (ex. --n), v:arguement (ex. 5)
    if len(all_args) == 11:
        arg_dict = {switch.lstrip("-"): argument for switch, argument in parse_args(args=all_args)}
    else:          # 6 args formatted like this: python ChiPSeqhelper.py --n=1 --narrowfile=narrowFH.txt ...(etc.)
        for arg in all_args[1:]:
            stripped_args_list = arg.lstrip("-").split("=")
            arg_dict[stripped_args_list[0]] = stripped_args_list[1]

    return arg_dict


def benchmarker(any_function):
    """
    :param any_function: (function to be timed)
    :return: String giving timing info for a given function
    TDLR: Decorator that takes a function as a parameter and outputs its timing info (benchmarking)
    """
    def timing_wrapper(*args, **kwargs):
        t1 = time.time()
        any_function(*args, **kwargs)
        t2 = time.time()
        return "Time it took to run {} function:{}\n".format(some_function.__name__, str(t2 - t1))
    return timing_wrapper


class ChipSeqPipeline(object):
    """
    Takes dictionary of validated arguments from check_args, sorts & selects MAC CLT output, and then calls IGV CLT
    """
    def __init__(self, args_dict):
        # Attributes
        self.args_dict = args_dict
        self.n = args_dict['n']
        self.genome = args_dict['genome']
        self.narrow_peaks = args_dict['narrowPeak_file']
        self.sortedNnarrowpeaks = self.narrow_peaks.split(".")[0] + "_SORTEDqValue_TOP50.narrowPeaks"  # Output Name
        self.treatment_bw = args_dict['treatmentBW_file']
        self.input_bw = args_dict['inputBW_file']
        self.output_folder = args_dict['output_folder']
        # Methods
        #self.validate()  # insert check args into this section
        self.run()

    def validate(self):   # insert check_args code into here
        """
        takes a list of __init__ attributes
        :return: boolean
        """
        pass

    def createIGVscript(self, inputBWfile, treatBWfile, sortedNarrowPeaksfile, genome, maxPanelheight=500,
                        padding=500, snapdirectory="snapstest"):
        """
        :param inputBWfile:      (ex. input.bw)
        :param treatBWfile:      (ex. treatment.bw)
        :param genome:        (ex. hg19, mm10, etc.)
        :param sortedNarrowPeaksfile:   # Used to grab chr# and start/stop positions
        :param maxPanelheight:  (default: 500px) # maximum screen shot panel size
        :param padding: (default: -/+ 500bp)     # left and right padding adjustments added to start/stop positions
        :param snapdirectory: (default:IGVsnaps) # Name of directory where the files will be stored
        ------------------------------------

        Example IGV scrpt:
        new
        snapshotDirectory IGV_Snapshots
        load test_alignments.bam
        genome hg19
        maxPanelHeight 500
        goto chr1:713167-714758
        snapshot chr1_713167_714758_h500.png
        goto chr1:713500-714900
        snapshot chr1_713500_714900_h500.png
        exit
        :return: (igv_batch_script.txt)
        """

        def get_workingDir():
            """
            Used for finding absolute paths
            :return: returns the current working directory
            """
            return os.getcwd()

        working_directory = get_workingDir()

        narrowFH = open(sortedNarrowPeaksfile, "r")
        outFH = open("igv_batch_script.txt", "w")

        # Writing out to the batch script file
        outFH.write("new\n")
        outFH.write("snapshotDirectory {}\n".format(snapdirectory))
        outFH.write("load {}\n".format(os.path.join(working_directory, treatBWfile)))
        outFH.write("load {}\n".format(os.path.join(working_directory, inputBWfile)))
        outFH.write("genome {}\n".format(genome))
        outFH.write("maxPanelHeight {}\n".format(maxPanelheight))
        for line in narrowFH:
            linelist = line.split()
            chr = linelist[0]  # chromosome number
            start = int(linelist[1]) - int(padding)  # start_position + some left padding
            stop = int(linelist[2]) + int(padding)  # stop_postion + some right padding
            outFH.write("goto {}:{}-{}\n".format(chr, start, stop))
            outFH.write("snapshot {}_{}_h{}.png\n".format(chr, start, stop, maxPanelheight))
        outFH.write("exit")
        outFH.close()
        narrowFH.close()


    def __run(self, command_list, pipe):
        """
        Private method used to run commands in shell
        When running pipe="yes", it runs:
        sort -k9nr,9 CHIP_Thpok_Biotin_vs_Input_Thpok_peaks.narrowPeak | head -50 > outputfile.narrowfile
        :param command_list:
        :pipe specify whether you want to pipe commands
        :return:
        """
        if pipe == "yes":
            p1 = sp.Popen(command_list, stdout=sp.PIPE)
            p2 = sp.Popen("head -{}".format(self.n).split(), stdin=p1.stdout, stdout=sp.PIPE)
            p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
            output = p2.communicate()[0]

            #fh = open("CHIP_Thpok_Biotin_vs_Input_Thpok_peaks_SORTEDqValue_TOP50.narrowPeak", "w")
            fh = open(self.sortedNnarrowpeaks, "w")
            fh.write(output)
            fh.close()

        elif pipe == "no":
            sp.Popen(command_list).wait()

    def run(self):
        """
        Sorts the Output of MACS (by decreasing q-value) & selects tops 'n' results
        calls IGV CLT, loops through results -> saves screenshots to an Output folder
        :return: Output folder with IGv N-peaks results
        """
        self.validate()
        self.__run("mkdir {}".format(self.output_folder).split(), "no")
        self.__run("sort -k9nr,9 CHIP_Thpok_Biotin_vs_Input_Thpok_peaks.narrowPeak".split(), "yes")
        self.__run("echo module load igv".split(), "no")
        self.createIGVscript(self.input_bw, self.treatment_bw, self.sortedNnarrowpeaks, genome=self.genome,
                             maxPanelheight=500, padding=500, snapdirectory=self.output_folder)

        #RUN THIS AFTER inserting testing.py: igv - m 20g - b igv_batch_script.txt
        #self.__run("echo Clean up the directory as needed-- rm any un-needed files!".split())

    def __str__(self):
        return "Parameters: {}".format(self.args_dict)



def main():
    """
    Pseudo-main method
    :return:
    """
    # Checking the Arguments pass through CL
    arg_list = sys.argv
    print(arg_list)
    args_dict = check_args(all_args=arg_list)
    print(args_dict)

    # Start working with interfacing into the Pipeline
    #for command in run_shell_commands(parsed_args=args_dict):
    #    sp.Popen("echo {}".format(command).split()).wait()  # Popen takes a list as parameter

    useChIPSeq = ChipSeqPipeline(args_dict)
    print(useChIPSeq)


if __name__ == "__main__":
    main()
