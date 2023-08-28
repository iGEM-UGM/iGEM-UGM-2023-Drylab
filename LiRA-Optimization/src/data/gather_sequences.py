import pandas as pd
import time
from selenium_profiles.webdriver import Chrome
from selenium_profiles.profiles import profiles
from selenium.webdriver.common.by import By  # locate elements
from selenium_profiles.utils.colab_utils import display, showscreen, show_html # virtual display

from webdriver_manager.chrome import ChromeDriverManager


def scrap_data_design(link, trial, driver):
  driver.get(f'{link}/{trial}')  # test client hints
  # Add a time delay to allow for page loading
  time.sleep(2)  # Adjust the delay as needed
  # dropdown sequences
  button = driver.find_elements(By.CSS_SELECTOR, "div.title")[1]
  button.click()
#   showscreen(driver) #for show driver overview

  #scrap the table of result
  tables = driver.find_elements(By.CLASS_NAME, "ui.table")
  table_data = []
  for tab in tables:
    rows = tab.find_elements(By.TAG_NAME, "tr")
    for row in rows:
        cells = row.find_elements(By.TAG_NAME, "td")
        row_data = [cell.text for cell in cells]
        table_data.append(row_data)

  df = pd.DataFrame(table_data, columns=["Domain_Strand", "Sequence"])
  return(df)

def get_design_data(link, driver):
  driver.get(link+'/0')
  # Find the menu element by its class name
  menu = driver.find_element(By.CLASS_NAME, "dropdown.icon")
  # Click the menu element
  menu.click()
  menu = driver.find_element(By.CSS_SELECTOR, "div.menu.transition.visible")
  # Find the items within the menu
  items = menu.find_elements(By.CLASS_NAME, "item")

  # Get the count of items
  item_count = len(items)
  print("Number of items:", item_count)

  # Create an empty DataFrame
  datas = pd.DataFrame()

  # Iterate over the items
  for i in range(item_count):
      # Scraping data using the scrap_data_design function
      data = scrap_data_design(link, i)

      n = len(data)  # The number of times to repeat the value
      result = [i+1] * n
      data['Trial'] = result
      # Append data to the datas DataFrame
      datas = datas.append(data, ignore_index=True)
  return(datas)

def main():
    url, email_address, password = input("Enter the URL: "), input("Enter the mail: "), input("Enter the password: ")
  
    chromedriver_path = ChromeDriverManager(version="90.0.4430.24").install()

    profile = profiles.Windows() # or .Android
    profile["cdp"]["cores"] = None # Chrome 90 doesn't allow emulating cores :(driver = mydriver.start(profile, uc_driver=False, executable_path=chromedriver_path)

    mydriver = Chrome(profile, executable_path=chromedriver_path)

    display = display()
    display.start_display()

    driver = mydriver.start()
    driver = mydriver.start()
    # Login
    driver.get('https://www.nupack.org')

    # Find the button by its CSS selector
    button = driver.find_element(By.CSS_SELECTOR, "a.item[href='/auth/log-in']")

    # Click the button
    button.click()

    # Find the input field by its class name
    input_field = driver.find_element(By.CLASS_NAME, "ui input")

    # Clear any existing text in the input field (optional)
    input_field.clear()

    # Enter the desired email address
    input_field.send_keys(email_address)
    input_element = driver.find_element(By.CLASS_NAME, "ui.fluid.icon.input")

    # Find the input field within the parent element
    password_field = input_element.find_element(By.TAG_NAME, "input")
    password_field.send_keys(password) #add your password here

    eye_icon = driver.find_element(By.CSS_SELECTOR, ".eye.fitted.link.icon")
    eye_icon.click()
    # Find the button
    button = driver.find_element(By.CLASS_NAME, "ui.button.orange")
    button.click()

    # start scrapping
    output_sequence = get_design_data(url)

    # Close the browser
    driver.quit()

    # save link into csv data downloaded automatically
    code = url.split('/')[-1]
    name = f"table_data_{code}.csv"
    output_sequence.to_csv('./'+name, index=False)
    
if __name__ == "__main__":
  main()
