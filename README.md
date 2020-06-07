# spectral_jaccard_similarity

## Setting up the environment 

The set of packages needed to run the experiments is provided [here](requirements.txt)

If one is using the [Anaconda package manager](https://www.anaconda.com/), one can set up the environment running the following commands

```
conda create -n sjs python=3.7 numpy=1.18.1
conda activate sjs
pip install -r requirements.txt
```

Then all the code in the repository can be run in the `sjs` environment created.

## Running code


To create filtered datasets run Bacterial_pipeline_part1.py and Bacterial_pipeline_part2.py. The `sjs` environment has all the dependencies needed to run them.

To use our pre-filtered datasets, in the `sjs` environment run,

```
python pipeline_wrapper.py --num_jobs <number-of-jobs> --datasets NCTC_ds.txt
```
