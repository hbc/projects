"""Retrieve top level 23andMe analyses for PGen participants into a spreadsheet.
"""
import argparse
import collections
import csv
import os

import requests
import flask
from flask import request

PORT = 5000
API_SERVER = "api.23andme.com"
BASE_API_URL = "https://%s/" % API_SERVER
BASE_CLIENT_URL = 'http://localhost:%s/' % PORT
REDIRECT_URI = '%sreceive_code/' % BASE_CLIENT_URL
DEFAULT_SCOPE = "basic analyses"
CLIENT_ID = os.environ["CLIENT_23ANDME_ID"]
CLIENT_SECRET = os.environ["CLIENT_23ANDME_SECRET"]

app = flask.Flask(__name__)

def get_user_info(headers):
    response = requests.get("%s1/user" % (BASE_API_URL),
                            headers=headers,
                            verify=False)
    res = response.json()
    return {"id": res["id"],
            "profiles": [x["id"] for x in res["profiles"] if x["genotyped"] is True]}

def get_pgen_patient_id(headers):
    user_info = get_user_info(headers)
    if len(user_info["profiles"]) == 1:
        profile_id = user_info["profiles"][0]
        response = requests.get("%s1/phenotypes/%s?phenotypes=bd_pgen_patient_id" %
                                (BASE_API_URL, profile_id),
                                headers=headers,
                                verify=False)
        res = response.json()
        return res["bd_pgen_patient_id"]
    else:
        print "Incorrect profile information", user_info
        return None

def summarize_analyses(access_token, pgen=False):
    """Provide a summary of available analysis traits for the given access token.
    """
    headers = {'Authorization': 'Bearer %s' % access_token}
    if pgen:
        profile_id = get_pgen_patient_id(headers)
        if profile_id is None:
            return
    else:
        profile_id = None
    Analysis = collections.namedtuple("Analysis", ["name", "attrs"])
    analyses = [Analysis("risks", ["risk", "population_risk"]),
                Analysis("carriers", ["mutations"]),
                Analysis("drug_responses", ["status"]),
                Analysis("traits", ["trait", "possible_traits"])]
    for analysis in analyses:
        response = requests.get("%s1/%s" % (BASE_API_URL, analysis.name),
                                headers=headers,
                                verify=False)
        if response.status_code == 200:
            res = response.json()
            if profile_id is None:
                assert len(res) == 1, "Need to handle users with multiple profiles"
                profile_id = res[0]["id"]
            for x in res[0][analysis.name]:
                info = [profile_id] + [x[attr] for attr in ["description"] + analysis.attrs]
                if len(analysis.attrs) < 2:
                    info.append("")
                if isinstance(info[-1], (list, tuple)):
                    info[-1] = ";".join(info[-1])
                yield info
        else:
            response.raise_for_status()

def write_summary_analyses(access_token, out_file):
    with open(out_file, "w") as out_handle:
        writer = csv.writer(out_handle)
        writer.writerow(["profile", "analysis", "result", "population"])
        for xs in summarize_analyses(access_token):
            writer.writerow([unicode(x).encode('ascii', errors='replace') for x in xs])

def write_pgen_analyses(access_tokens, out_file, starttoken=None):
    problem_file = "%s-problems.txt" % os.path.splitext(out_file)[0]
    problem_handle = open(problem_file, "w")
    with open(out_file, "w") as out_handle:
        writer = csv.writer(out_handle)
        writer.writerow(["patientid", "analysis", "result", "population"])
        print "Processing ", len(access_tokens), " tokens"
        start = False
        for i, access_token in enumerate(access_tokens):
            print i, access_token
            if starttoken is None or starttoken == access_token:
                start = True
            if start:
                analyses = list(summarize_analyses(access_token, pgen=True))
                if len(analyses) == 0:
                    print "No genotype data"
                    problem_handle.write("%s %s\n" % (i + 1, access_token))
                else:
                    for xs in summarize_analyses(access_token, pgen=True):
                        writer.writerow([unicode(x).encode('ascii', errors='replace') for x in xs])

def get_pgen_access_tokens():
    """Retrieve PGen access tokens from specialized endpoint at 23andme.
    """
    parameters = {
        'client_id': CLIENT_ID,
        'client_secret': CLIENT_SECRET}
    response = requests.post(
        "%s%s" % (BASE_API_URL, "access_tokens/"),
        data=parameters,
        verify=False
    )
    if response.status_code == 200:
        tokens_json = response.json()
        return tokens_json["access_tokens"]
    else:
        response.raise_for_status()

@app.route('/')
def index():
    auth_url = "%sauthorize/?response_type=code&redirect_uri=%s&client_id=%s&scope=%s" % (
        BASE_API_URL, REDIRECT_URI, CLIENT_ID, DEFAULT_SCOPE)
    return flask.render_template('index.html', auth_url=auth_url)

@app.route('/receive_code/')
def receive_code(code=None):
    parameters = {
        'client_id': CLIENT_ID,
        'client_secret': CLIENT_SECRET,
        'grant_type': 'authorization_code',
        'code': code or request.args.get('code'),
        'redirect_uri': REDIRECT_URI,
        'scope': DEFAULT_SCOPE,
    }
    response = requests.post(
        "%s%s" % (BASE_API_URL, "token/"),
        data = parameters,
        verify=False
    )
    print response.text
    if response.status_code == 200:
        access_json = response.json()
        print "Access token", access_json["access_token"]
        return access_json["access_token"]
    else:
        response.raise_for_status()

if __name__ == '__main__':
    parser = argparse.ArgumentParser("Batch retrieve or get access tokens from 23andMe")
    parser.add_argument("--token", help="Single user access token to retrieve")
    parser.add_argument("--out", help="Output file to write CSV results")
    parser.add_argument("--pgen", action='store_true', default=False,
                        help="Run retrieval of PGen specific results")
    parser.add_argument("--starttoken", help="Token to start pgen analysis at")
    args = parser.parse_args()
    if args.pgen:
        tokens = get_pgen_access_tokens()
        write_pgen_analyses(tokens, args.out, args.starttoken)
    elif args.token is not None:
        print "Retrieving results for", args.token
        write_summary_analyses(args.token, args.out)
    else:
        print "A local client for the Personal Genome API is now initialized."
        app.run(debug=False, port=PORT)
