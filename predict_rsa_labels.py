### DEPENDENCIES

import os
import pickle
import numpy as np
import pandas as pd
from tensorflow import keras
import argparse

### FUNCTIONS

###### AUXILIARY FUNCTIONS

def load_pickle(f_in):
    """
    loads pickle
    """
    with open(f_in, "rb") as f:
        data = pickle.load(f)
    return data

def get_RSA_vectors(rsa_profs_filt):
    """
    Given a ddictionary of RSA profiles (arrays),
    generates a dataframe of RSA vectors representing
    a binding site.
    """
    max_len = max([len(v) for v in rsa_profs_filt.values()])
    bs_vectors = []
    bs_vectors_dict = {}
    for bs_id, rsa_sig in rsa_profs_filt.items():
        rsa_sig_len = len(rsa_sig)
        rsa_range_prop = [0 for i in range(10)]  # now let us change to 10
        for rsa in rsa_sig:
            prop_i = int(rsa / 10)  # 10 RSA BINS: b1 = [0,10), b2 = [10, 20), ... b10 = [90, MAX)
            if prop_i > 9:  # if greater than 100, put in 10th bin
                prop_i = 9
            rsa_range_prop[prop_i] += 1
        rsa_range_prop = [round(i / rsa_sig_len, 3) for i in rsa_range_prop]
        rsa_range_prop.insert(0, rsa_sig_len / max_len)  # ADDING BINDING SITE SIZE RELATIVE TO MAX SITE SIZE (IN THIS CASE 40)
        bs_vectors.append(rsa_range_prop)
        bs_vectors_dict[bs_id] = rsa_range_prop

    vector_df = pd.DataFrame(bs_vectors, index=list(rsa_profs_filt.keys()))  # obtaining RSA vectors, which are the 11-element features used for the machine learning

    return vector_df

###### MAIN FUNCTION

def main(input_dir):
    """
    This function will add RSA Cluster labels and RSA-based functional scores
    to a summary binding site table dataframe.
    """
    results_dir = f'./OUT/{input_dir}/results'
    bss_data_out = os.path.join(results_dir, f'{input_dir}_bss_table.pkl')
    bss_data_out_RSA = os.path.join(results_dir, f'{input_dir}_bss_RSA_table.pkl')
    rsa_profs_out = os.path.join(results_dir, f'{input_dir}_bss_RSA_profiles.pkl')

    model_path = "./OTHER/RSA_pred_model.h5"

    # Load data
    bss_data = pd.read_pickle(bss_data_out)
    rsa_profs = load_pickle(rsa_profs_out)
    rsa_profs_filt = {k: v for k, v in rsa_profs.items() if len(v) > 0}

    # Get RSA vectors
    vector_df = get_RSA_vectors(rsa_profs_filt)

    # Load the model
    final_model = keras.models.load_model(model_path)

    # Predict
    final_preds = final_model.predict(x=vector_df, batch_size=27, verbose=0)
    rounded_predictions = np.argmax(final_preds, axis=-1)
    site_names = vector_df.index.tolist()
    
    clab_dict = {}
    for i in range(len(vector_df)):
        clab_dict[site_names[i]] = rounded_predictions[i]+1

    # Calculate functional scores
    func_scores_dict = {}
    for i in range(len(final_preds)):
        func_score = ((0.52 * final_preds[i][0]) + (0.18 * final_preds[i][1]) + (0.05 * final_preds[i][2]) + (0.04 * final_preds[i][3]))
        func_scores_dict[site_names[i]] = round(func_score, 2)

    # Add data to DataFrame
    bss_data["Cluster"] = bss_data.ID.map(clab_dict)
    bss_data["FS"] = bss_data.ID.map(func_scores_dict)

    # Save the modified DataFrame
    bss_data.to_pickle(bss_data_out_RSA)
    
### RUNNING SCRIPT

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This script predicts RSA cluster labels and calculates RSA-based functional score (FS)')
    parser.add_argument('input_dir', type=str, help='This is the Input ID or Job ID, i.e., name of the directory where the binding site table resides.')
    args = parser.parse_args()

    main(args.input_dir)

