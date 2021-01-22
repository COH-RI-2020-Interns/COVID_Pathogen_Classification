from selenium import webdriver
from os import getcwd
from selenium.webdriver.common.keys import Keys
from bs4 import BeautifulSoup
import re
import pandas as pd

getcwd()
DRIVER_PATH = getcwd() + '/chromedriver'
driver = webdriver.Chrome(executable_path=DRIVER_PATH)
driver.get('https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=2559587')
family = driver.find_element_by_xpath("//input")
