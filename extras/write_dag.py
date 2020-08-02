import argparse

def uber_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-S", "--stacking", action="store_true", default=False,
                        help="Use this option to create stacking dag")

    return parser

def parser_compute_evidence():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dagfile", action="store", help="Name of the output DAG file")
    parser.add_argument("-i", "--input", action="store", help="Name of posterior file")
    parser.add_argument("-p", "--prior", action="store", help="Name of the prior file")
    parser.add_argument("-o", "--output", action="store", help="Prefix of the output file names")
    parser.add_argument("-n", "--num", action="store", type=int, help="Number of jobs")
    parser.add_argument("-t", "--target", action="store", help="Name of target EoS")
    parser.add_argument("-r", "--reference", action="store", help="Name of reference EoS")
    parser.add_argument("-T", "--trials", action="store", type=int, help="Number of trials")
    parser.add_argument("-j", "--joinsubfilename", action="store", help="Name of json joining subfile")
    parser.add_argument("-s", "--subfilename", action="store", help="Name of bayes-factor computation subfile")
    return parser

def parser_stacking():
    parser = argparse.ArgumentParser(parents=[random_parser])
    parser.add_argument("-d", "--dagfile", action="store", help="Name of the output DAG file")
    parser.add_argument("--posteriorlist", action="store",
                        help="A text file with list of posterior samles for stacking")
    parser.add_argument("--priorlist", action="store",
                        help="A text file with List of priors samples for stacking")
    parser.add_argument("-o", "--output", action="store", help="Prefix of the output file names")
    parser.add_argument("-n", "--num", action="store", type=int, help="Number of jobs")
    parser.add_argument("-t", "--target", action="store", help="Name of target EoS")
    parser.add_argument("-r", "--reference", action="store", help="Name of reference EoS")
    parser.add_argument("-T", "--trials", action="store", type=int, help="Number of trials")
    parser.add_argument("-S", "--stacksubfilename", action="store",
                        help="Name of stacking computation subfile")
    return parser


def main(args=None):
    import os

    P = uber_parser()
    uber_args = P.parse_args()

    ## Create directories ##
    os.system("mkdir -p output")
    os.system("mkdir -p error")
    os.system("mkdir -p logs")
    os.system("mkdir -p {}".format(args.output))

    if uber_args.stacking:
        p = parser_stacking()
        args = p.parse_args()
        stack_subfile_txt = '''universe = vanilla
executable = ./comp_stacked_evidence
arguments = "--input {} --prior {} --target {} --reference {} --nums 40 --output {}"
output = output_$(macrotag).stdout
error = error_$(macrotag).stderr
log = logfile_$(macrotag).log
getenv = True
accounting_group = ligo.prod.o2.cbc.pe.lalinferencerapid
queue 1
'''.format(posterior_files_list, prior_files_list.dat, args.target, args.reference, args.output)

        with open(args.stacksubfilename, 'w') as f:
            f.writelines(stack_subfile_txt)

    else:
        p = parser_compute_evidence()
        args = p.parse_args()

        ### Write comp_eos_evidence sub file ###
        subfile_txt = '''universe = vanilla
executable = ./comp_eos_evidence
arguments = "--input $(macroinput) --prior $(macroprior) --target {} --reference {} --nums {} --output $(macrooutput)"
output = output/output_$(macrotag).stdout
error = error/error_$(macrotag).stderr
log = logs/logfile_$(macrotag).log
getenv = True
accounting_group = ligo.prod.o2.cbc.pe.lalinferencerapid
queue 1
'''.format(args.target, args.reference, args.trials)

        with open(args.subfilename, 'w') as f:
            f.writelines(subfile_txt)

        ### Write join jsons sub file ###
        subfile_txt ='''universe = vanilla
executable = ./combine_results
arguments = "--jsons {} --target {} --reference {} --output {}.json"
output = output/output_join.stdout
error = error/error_join.stderr
log = logs/logfile_join.log
getenv = True
accounting_group = ligo.prod.o2.cbc.pe.lalinferencerapid
queue 1
'''.format(args.output, args.target, args.reference, args.output)

        with open(args.joinsubfilename, 'w') as f:
            f.writelines(subfile_txt)

        ### Write the DAG file ###
        with open(args.dagfile, 'w') as f:
            text_parent_child = "PARENT"
            for ii in range(args.num):
                text_parent_child += " " + str(ii)
                text = '''JOB {} {}\nVARS {} macroinput="{}"\nVARS {} macroprior="{}"\nVARS {} macrooutput="{}/{}_{}.json"\nVARS {} macrotag="{}"\n\n'''.format(ii, args.subfilename, ii, args.input, ii, args.prior, ii, args.output, args.output, ii, ii, ii)
                f.writelines(text)

            text_join = "JOB {} {}\n".format(ii+1, args.joinsubfilename)
            f.writelines(text_join)
            text_parent_child += " CHILD " + str(ii+1)
            f.writelines(text_parent_child)

if __name__ == "__main__":
    main()

### Write stacking sub file ###
#if args.posteriors:
#    stack_subfile_txt = '''universe = vanilla
#executable = ./comp_stacked_evidence
#arguments = "--input $(macroposteriors) --prior (macropriorfiles) --labels $(macrolabels) --target {} --reference {} --nums {} --output $(macrooutput)"
#output = output_$(macrotag).stdout
#error = error_$(macrotag).stderr
#log = logfile_$(macrotag).log
#getenv = True
#accounting_group = ligo.prod.o2.cbc.pe.lalinferencerapid
#queue 1
#'''.format(args.target, args.reference, args.trials)

#    with open(args.stacksubfilename, 'w') as f:
#        f.writelines(stack_subfile_txt)

