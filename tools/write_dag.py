import argparse


def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dagfile", action="store",
                        help="Name of the output DAG file")
    parser.add_argument("-i", "--input", action="store", default=None,
                        help="Name of posterior file")
    parser.add_argument("--posteriorlist", action="store", default=None,
                        help="Text file with list of posteriors for stacking")
    parser.add_argument("-p", "--prior", action="store", default=None,
                        help="Name of the prior file")
    parser.add_argument("--priorlist", action="store", default=None,
                        help="Text file with List of priors for stacking")
    parser.add_argument("-o", "--output", action="store",
                        help="Prefix of the output file names")
    parser.add_argument("-n", "--num", action="store", type=int,
                        help="Number of jobs")
    parser.add_argument("-t", "--target", action="store",
                        help="Name of target EoS")
    parser.add_argument("-r", "--reference", action="store",
                        help="Name of reference EoS")
    parser.add_argument("-T", "--trials", action="store", type=int,
                        help="Number of trials")
    parser.add_argument("-j", "--joinsubfilename", action="store",
                        help="Name of json joining subfile")
    parser.add_argument("-s", "--subfilename", action="store",
                        help="Name of bayes-factor computation subfile")
    parser.add_argument("-S", "--stacksubfilename", action="store",
                        help="Name of stacking computation subfile")
    return parser


def main(args=None):
    import os
    import sys

    p = parser()
    args = p.parse_args()

    ## Create directories ##
    output_dir = "output_{}".format(args.output)
    os.system("mkdir -p {}".format(output_dir))
    error_dir = "error_{}".format(args.output)
    os.system("mkdir -p {}".format(error_dir))
    logs_dir = "logs_{}".format(args.output)
    os.system("mkdir -p {}".format(logs_dir))
    os.system("mkdir -p {}".format(args.output))

    ## Sanity checks ##
    # if (args.input is None) ^ (args.posteriorlist is None):
    #     print('Either --input or --posteriorlist arguments should be provided')
    #     sys.exit(0)


    if args.posteriorlist:
        stack_subfile_txt = '''universe = vanilla
executable = ./comp_stacked_evidence
arguments = "--input {} --prior {} --target {} --reference {} --nums {} --output $(macrooutput)"
output = {}/output_$(macrotag).stdout
error = {}/error_$(macrotag).stderr
log = {}/logfile_$(macrotag).log
getenv = True
accounting_group = ligo.prod.o2.cbc.pe.lalinferencerapid
queue 1
'''.format(args.posteriorlist, args.priorlist, args.target, args.reference, args.trials, output_dir, error_dir, logs_dir)

        with open(args.stacksubfilename, 'w') as f:
            f.writelines(stack_subfile_txt)

        join_subfile_txt = '''universe = vanilla
executable = ./combine_stacked_results
arguments = "--json {} --output {}"
output = {}/output_$(macrotag).stdout
error = {}error_$(macrotag).stderr
log = {}/logfile_$(macrotag).log
getenv = True
accounting_group = ligo.prod.o2.cbc.pe.lalinferencerapid
queue 1
'''.format(args.output, args.output, output_dir, error_dir, logs_dir)

        with open(args.joinsubfilename, 'w') as f:
            f.writelines(join_subfile_txt)

        ### Write the DAG file ###
        with open(args.dagfile, 'w') as f:
            text_parent_child = "PARENT"
            for ii in range(args.num):
                text_parent_child += " " + str(ii)
                text = '''JOB {} {}\nVARS {} macrooutput="{}/{}_{}.json"\nVARS {} macrotag="{}"\n\n'''.format(ii, args.stacksubfilename, ii, args.output, args.output, ii, ii, ii)
                f.writelines(text)

            text_join = '''JOB {} {}\nVARS {} macrotag="{}"\n'''.format(ii+1, args.joinsubfilename, ii+1, args.output)
            f.writelines(text_join)
            text_parent_child += " CHILD " + str(ii+1)
            f.writelines(text_parent_child)


    if args.input:
        ### Write comp_eos_evidence sub file ###
        subfile_txt = '''universe = vanilla
executable = ./comp_eos_evidence
arguments = "--input $(macroinput) --prior $(macroprior) --target {} --reference {} --nums {} --output $(macrooutput)"
output = {}/output_$(macrotag).stdout
error = {}/error_$(macrotag).stderr
log = {}/logfile_$(macrotag).log
getenv = True
accounting_group = ligo.prod.o2.cbc.pe.lalinferencerapid
queue 1
'''.format(args.target, args.reference, args.trials, output_dir, error_dir, logs_dir)

        with open(args.subfilename, 'w') as f:
            f.writelines(subfile_txt)

        ### Write join jsons sub file ###
        subfile_txt ='''universe = vanilla
executable = ./combine_results
arguments = "--jsons {} --target {} --reference {} --output {}.json"
output = {}/output_join.stdout
error = {}/error_join.stderr
log = {}/logfile_join.log
getenv = True
accounting_group = ligo.prod.o2.cbc.pe.lalinferencerapid
queue 1
'''.format(args.output, args.target, args.reference, args.output, output_dir, error_dir, logs_dir)

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
