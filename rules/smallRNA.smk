# function to access logs more easily
def return_log_smallrna(sample_name, step, paired):
    return os.path.join(REPO_FOLDER,"smallRNA","logs",f"tmp__{sample_name}__{step}__{paired}.log")

