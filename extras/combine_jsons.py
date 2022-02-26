import os
import json
import numpy as np
import argparse

def join_json_files(list_of_jsons, nametag="model"):
    '''
    This helper function joins the JSON files from the output of the Bayes
    factor computation methods in Model_selection and Stacking into one single
    JSON file that can be directly used for further analysis.

    list_of_jsons :: A list of all the JSON files that needs to be combined.
    nametag :: Use this option to create a name-tag for the resulting JSON
    files
    '''
    alldata = []
    reference_models = []
    for json_file in list_of_jsons:
        if os.path.exists(json_file):
            with open(json_file, 'r') as f:
                try:
                    thisdata = json.load(f)
                    reference_models.append(thisdata['ref_eos'])
                    alldata.append(thisdata)
                except TypeError:
                    print("{} is not formatted correctly".format(json_file))
                    quit()
                except KeyError:
                    print("{} is not formatted correctly".format(json_file))
                    quit()
        else:
            print("Could not find {}".format(json_file))

    # Keep unique names of the equation of state models
    reference_models = np.unique(np.array(reference_models)).tolist()

    for model in reference_models:
        combined_dict = {}
        for data in alldata:
            if model == data['ref_eos']:
                value = {'joint_bf': data['joint_bf'],
                         'joint_bf_array': data['joint_bf_array'],
                         'all_bf': data['all_bf'],
                         'all_bf_err': data['all_bf_err']}
                combined_dict[data['target_eos']] = value

        filename = 'bayes_factors_against_' + model + '_' + nametag + '.json'
        with open(filename, 'w') as f:
            json.dump(combined_dict, f, indent=2, sort_keys=True)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-C", "--cache", action="store",
                        help="A cache of JSON file names that needs to be joined")
    parser.add_argument("-t", "--tag", action="store",
                        help="name-tag of the output JSON file")

    args = parser.parse_args()

    file_names = list(np.loadtxt(args.cache, dtype=str))

    join_json_files(file_names, nametag=args.tag)
