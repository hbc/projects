## PGen 23andme response project

Retrieve 23andme analysis results for PGen project. 

### Setup

    cd pgen_23andme
    virtualenv .
    source bin/activate
    pip install -r requirements.txt
    
## Usage

Before working, activate the local python with installations:

    source bin/activate
   
To run a small web server to get your own personal access token:
    
    python pgen_retrieval.py
    google-chrome http://127.0.0.1:5000/
    
To retrieve results for a single user access token:
    
    python pgen_retrieval.py --token=example_token --out=output.csv

Framework code taken from 23andMe's flask example:

https://github.com/23andMe/api-example-flask
