# spectral_jaccard_similarity

# Setting up the environment 

The set of packages needed to run the experiments is provided [here](requirements.txt)

If one is using the [Anaconda package manager](https://www.anaconda.com/), one can set up the environment running the following commands

```
conda create -n sjs python=3.7
conda activate sjs
pip install -r requirements.txt
```

Then all the code in the repository can be run in the `sjs` environment created.


## Dependencies
- mmh3 hash library
- SimpleFastaParser from Bio.SeqIO.FastaIO
- tqdm progress bars

## Instructions
To create filtered datasets run Bacterial_pipeline_part1.py and Bacterial_pipeline_part2.py.

To use our pre-filtered datasets first run makeFolders.py then run pipeline_wrapper.py. All the code is encapsulated in there.
