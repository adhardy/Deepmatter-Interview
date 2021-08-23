import os
import pathlib

import xml.dom.pulldom as pulldom
import xml.etree.ElementTree as ET

import matplotlib.dates as mdates
import matplotlib.pyplot as plt

import pandas as pd
from typing import List

def check_file_extension(file_path: str, extension: str): # -> bool:
    """Checks if a file has a specific extension"""
    return pathlib.Path(file_path).suffix == extension

def parse_pccr(pccr_dir: List[str]):
    pccr_file_ext = ".pcrr"
    pccr_files = os.listdir(pccr_dir) # get all files in the directory
    recipe_information = []

    for pccr_file in pccr_files: # iterate over each file

        # === Check for correct file type ===
        if not check_file_extension(pccr_file, pccr_file_ext):
            print(f"{pccr_file} is not a .pccr file, skipping.")
            continue
        
        pccr_file_path = os.path.join(pccr_dir,pccr_file) # get the full path

        # === Parse XML ===
        doc = pulldom.parse(pccr_file_path) # these xml files are very large, use pulldom to extract the parts we need 
        for event, node in doc:
            if event == pulldom.START_ELEMENT and node.tagName == "recipe_version": # all required information is within this tag
                doc.expandNode(node) # expand the node so we can parse it with elementree and xpath
                root = ET.fromstring(node.toxml())  # load the xml into elementree

                # use xpath to find recipe name and author
                name = root.find("./recipe").attrib["name"]
                author = root.find("./pcml/meta/author").text

                # use xpath to find all operation tags, and then count them
                operations = len(root.findall("./pcml/step/group/operation"))
                
                recipe_information.append({"filename":pccr_file, "name": name, "author": author, "operations": operations})

                break # we don't need any more information, stop processing the file
    
    return recipe_information

def extract_sensor_data(pcrr_file_path: str,
                        sensor_name: str,
                        stop_after_first: bool = False): # -> pd.DataFrame: # adding the return type hints messes up the syntax highlightning in notebooks in VSCode

    """Extracts values and timestamps from the given pcrr file and sensor.
    pccr_file_path: Path to the pcrr file
    sensor_name: The name of the sensor to read data for
    stop_after_first: Defaults to False. Set True to stop parsing the document after the first instance of the sensor is found.
    """
    
    # === Check for correct file type ===
    pccr_file_ext = ".pcrr"
    if not check_file_extension(pcrr_file_path, pccr_file_ext):
        print(f"{pcrr_file_path} is not a .pccr file, aborting.")
        return None

    sensor_data_value = []
    sensor_data_timestamp = []

    # === XML Parsing ===
    doc = pulldom.parse(pcrr_file_path) # these xml files are very large, use pulldom to extract the parts we need
    for event, node in doc:
        if event == pulldom.START_ELEMENT:
            # could extend this here by iterating over multiple sensors, meaning we would only have to process the xml file once rather than n-sensor times
            if node.tagName == "sensor_data" and node.getAttribute("name") == sensor_name: # find the sensor_data tag with the name attribute of the sensor we're interested in
                doc.expandNode(node) # expand the node so we can parse it with elementree and xpath
                root = ET.fromstring(node.toxml()) # load the xml into elementree

                # find all instances of the sensor_data_record tag within this node, record the values and timestamps
                sensor_data_value += [x.text for x in root.findall("./sensor_data_record/value")] 
                sensor_data_timestamp += [x.text for x in root.findall("./sensor_data_record/timestamp")]

                if stop_after_first: # sensor data is present in mutliple nodes. If set, stop parsing after the first instance of the target tag, useful for debugging/demo purposes as it is otherwise quite slow because of the large file size
                    break

    # === Create Dataframe ===
    df = pd.DataFrame(data={sensor_name: sensor_data_value, "timestamp": sensor_data_timestamp}) #convert to a dataframe
    df = df.astype({sensor_name: float, "timestamp": "datetime64[ns]"}) # set datatypes
    df.set_index("timestamp", inplace=True) # set the timestamp as the index
    
    return df

def plot_sensor_data(df: pd.DataFrame, 
                    y_label: str = None, 
                    df2: pd.DataFrame = None, 
                    y2_label: str = None): # -> plt.Figure: # adding the return type hints messes up the syntax highlightning in notebooks in VSCode
    """Plots sensor data extracted from pcrr file by 'extract_sensor_data'
    df1: Dataframe to plot on first y axis
    y_label: Label for y axis
    df2: Dataframe to plot on second y axis
    y2_label: Label for second y axis"""
    
    fig, axs = plt.subplots()
        
    df.plot(ax=axs,
            xlabel="Date/Time", 
            ylabel=y_label) 
    

    if df2 is not None:
        ax2 = axs.twinx() # create a second y axis
        df2.plot(ax=ax2,
                ylabel=y2_label,
                color="r")
                
        # if using multiple axis, need to curate legend manually
        lines, labels = axs.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        axs.legend(lines + lines2, labels + labels2)
    else:
        # if not using multiple axis, just use the legend from the first axis
        axs.legend()

    timestamp_format = mdates.DateFormatter('%d/%m %H:%M') # set the datetime format of the x-axis 
    axs.xaxis.set_major_formatter(timestamp_format)    

    return fig