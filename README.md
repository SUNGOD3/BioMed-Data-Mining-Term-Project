# [BioMed Data Mining Term Project](https://github.com/SUNGOD3/BioMed-Data-Mining-Term-Project)
Implement a program to predict the secondary structure commonly shared in an RNA 
family.

## Method Overview

* Uses covariation analysis to identify conserved base pairs across the RNA family
* Considers both canonical (A-U, G-C) and non-canonical (G-U) base pairs
* Implements a scoring system based on the frequency of valid base pairs at each position


## Key Components

* Sequence reading and alignment handling
* Covariation matrix calculation
* Structure prediction based on covariation scores
* Output formatting in the required notation


## Algorithm Steps

* Read FASTA sequences
* Create/handle multiple sequence alignment
* Calculate covariation scores between all positions
* Identify likely base pairs based on covariation scores
* Generate dot-bracket notation for the predicted structure


## Performance Considerations

* The program tracks TP, TN, FP, and FN implicitly through the covariation scoring
* Uses a minimum score threshold (default 0.75) to control prediction stringency
* Implements basic pseudoknot avoidance

## How to Run

1. **Clone the repository**:
    ```sh
    git clone https://github.com/SUNGOD3/BioMed-Data-Mining-Term-Project.git
    cd BioMed-Data-Mining-Term-Project
    ```

2. **Install dependencies**:
    Make sure you have all the necessary dependencies installed. You can use the following command to install them:
    ```sh
    pip install -r requirements.txt
    ```

3. **Prepare your input data**:
    Ensure your input RNA sequences are in FASTA format and placed in the appropriate directory.
    e.g. RNA_Dataset1.txt

4. **Run the program & View the results**:
    ```sh
    python main.py "Dir of your input file"
    ```
    The execution results will be displayed on the terminal