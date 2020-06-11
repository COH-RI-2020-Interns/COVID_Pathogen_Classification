# Automation of Fasta Sequence Extraction on SourceForge

## Description:

This is a Python script developed to provide ease in the process of collecting .fasta files by Randhawa et. al. from SourceForge.

## Instructions:

### 1. Pull this script to your branch.

```
$ git checkout [your branch name]
$ git pull origin fasta_scraping
```

### 2. Install all the requirements for the Python script.

- If Python 3 is your default Python:

```
$ pip install -r src/requirements.txt
```
- If Python 3 is not your default Python:

```
$ pip3 install -r src/requirements.txt
```

### 3. Run the Script in your Terminal

- If Python 3 is your default Python:

```
$ python src/covid_sourceforge_scraper.py
```
- If Python 3 is not your default Python:

```
$ python3 src/covid_sourceforge_scraper.py
```

You should see the data download immediately in the data folder under the specified virus family input. 
