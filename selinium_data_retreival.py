from selenium import webdriver
from selenium.webdriver.common.keys import Keys
import time
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
from os import getcwd
driver=webdriver.Chrome(getcwd() + "/chromedriver")
driver.get("https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide")
paths="//*[@id='filter-toggle-VirusLineage_ss']"
select=driver.find_element_by_xpath(paths).click()
enter_virus=driver.find_element_by_css_selector("#typeahead-input-VirusLineage_ss").send_keys("riboviria")
time.sleep(1)
select_virus=driver.find_element_by_css_selector(".flyout-panel .usa-sidenav-list.typeahead-results button");
select_virus.click()
time.sleep(0.5)
select_complete=driver.find_element_by_xpath("//*[@id='Completeness_s-panel-toggle']").click()
time.sleep(1)
select_complete=driver.find_element_by_xpath("//*[@id='Completeness_s-filter-panel']/uswds-ncbi-app-facet-item-checkbox[1]/div/label").click()
select_target_date=driver.find_element_by_xpath("//*[@id='filter-toggle-CreateDate_dt']").click()
enter_start_date=driver.find_element_by_xpath("//*[@id='CreateDate_dt_From']").send_keys("1/1/2016")
time.sleep(0.3)
enter_end_date=driver.find_element_by_id("CreateDate_dt_To")
time.sleep(0.3)
enter_end_date.send_keys("null")
time.sleep(0.3)
enter_end_date.clear()
enter_end_date.send_keys("12/31/2018")
subit_date=driver.find_element_by_xpath("//*[@id='submit-CreateDate_dt']").click()
time.sleep(0.3)
i=0
try:
    for b in range(1,6):
        for a in range(1,201):
            time.sleep(5)
            check_box=driver.find_element_by_css_selector("#DataTables_Table_0 > tbody > tr:nth-child("+str(a)+") > td:nth-child(1)").click()
            time.sleep(5)
            download=driver.find_element_by_xpath("//*[@id='cmscontent']/section/uswds-ncbi-app-root/uswds-ncbi-app-report/div/div[2]/uswds-ncbi-app-tool-button/div/uswds-ncbi-app-download-module/button").click()
            time.sleep(0.3)
            download_next = driver.find_element_by_xpath("//*[@id='1']/div/button").click()
            time.sleep(0.3)
            download_next_next = driver.find_element_by_xpath("//*[@id='2']/div/span[1]/button").click()
            time.sleep(0.3)
            final_download = driver.find_element_by_xpath("//*[@id='3']/div/span[1]/button[2]").click()
            time.sleep(5)
            checkbox = WebDriverWait(driver, 120).until(EC.element_to_be_clickable((By.CSS_SELECTOR,"#DataTables_Table_0 > tbody > tr:nth-child("+str(a)+") > td:nth-child(1)" ))).click()

            time.sleep(2)

            i+=1
        if i==1000:
            break
        else:
            next_page=driver.find_element_by_xpath("//*[@id='cmscontent']/section/uswds-ncbi-app-root/uswds-ncbi-app-report/div/div[2]/uswds-ncbi-app-report-data/div/div[3]/div/div/uswds-ncbi-app-results-pagination/form/div/button[2]").click()
            continue
except:
    print("No more records")
