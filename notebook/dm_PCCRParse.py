import os
import pathlib

import xml.dom.pulldom as pulldom
import xml.etree.ElementTree as ET

import matplotlib.dates as mdates
import matplotlib.pyplot as plt

import pandas as pd
from typing import List,Any, Union

pccr_file_ext = ".pcrr"

def check_file_extension(file_path: str, extension: str): # -> bool:
    """Checks if a file has a specific extension"""
    if pathlib.Path(file_path).suffix != extension:
        raise TypeError(f"{file_path} is not a .pccr file.")

def none_or_equal(value1:any, value2: any):# -> bool:
    """Returns True if the values are equal or if either value is None. Otherwise returns False."""
    if value1 is None or value2 is None or value1 == value2:
        return True
    return False

def get_XML_root(pccr_file_path: str, tagname: str, attribute: str = None, attribute_value: Any = None):
    """Generator: yields XML trees matching the given tag and attribute."""

    doc = pulldom.parse(pccr_file_path) # these xml files are very large, use pulldom to extract the parts we need 
    for event, node in doc:
        if event == pulldom.START_ELEMENT and none_or_equal(node.tagName, tagname) and none_or_equal(node.getAttribute(attribute), attribute_value):
            doc.expandNode(node) # expand the node so we can parse it with elementree and xpath
            yield ET.fromstring(node.toxml())  # load the xml into elementree

def get_pccr_metadata(pccr_dir: str):
    """Extract name, author and number of operations from a pccr files located in a directory."""

    pccr_files = os.listdir(pccr_dir) # get all files in the directory
    recipe_information = []

    for pccr_file in pccr_files: # iterate over each file

        check_file_extension(pccr_file, pccr_file_ext) # Check for correct file type
        pccr_file_path = os.path.join(pccr_dir,pccr_file) # get the full path

        for root in get_XML_root(pccr_file_path, tagname="recipe_version"):
            
            name = root.find("./recipe").attrib["name"] # use xpath to find recipe name and author
            author = root.find("./pcml/meta/author").text
            operations = len(root.findall("./pcml/step/group/operation")) # use xpath to find all operation tags, and then count them
            recipe_information.append(RecipeInformation(name, author, operations, pccr_file)) # append the recipe information to a list
            break # we don't need any more information, stop processing the file
    
    return recipe_information

class RecipeInformation():
    """Class to store recipe information"""
    def __init__(self, name: str, author: str, operations:int, filename: str):
        self.name = name
        self.author = author
        self.operations = operations
        self.filename = filename

    def __repr__(self):
        return f"{self.filename}: {self.name} by {self.author} with {self.operations} operations."


def extract_sensor_data(pccr_file_path: str,
                        sensors: Union[str, List[str]],
                        stop_after_first: bool = False): # -> pd.DataFrame: # adding the return type hints messes up the syntax highlightning in notebooks in VSCode

    """Extracts values and timestamps from the given pcrr file and sensor.
    pccr_file_path: Path to the pcrr file
    sensor_name: The name of the sensor to read data for
    stop_after_first: Defaults to False. Set True to stop parsing the document after the first instance of the sensor is found.
    """

    if not isinstance(sensors, list):
        sensors = [sensors]

    check_file_extension(pccr_file_path, pccr_file_ext) # Check for correct file type 
    
    data = {} # dictionary to hold the extracted data
    # === XML Parsing ===
    for sensor_name in sensors:
        sensor_data_value = []
        sensor_data_timestamp = []
        for root in get_XML_root(pccr_file_path, tagname="sensor_data", attribute="name", attribute_value=sensor_name):

            # find all instances of the sensor_data_record tag within this node, record the values and timestamps
            sensor_data_value += [x.text for x in root.findall("./sensor_data_record/value")] 
            sensor_data_timestamp += [x.text for x in root.findall("./sensor_data_record/timestamp")]

            if stop_after_first: # sensor data is present in mutliple nodes. If set, stop parsing after the first instance of the target tag, useful for debugging/demo purposes as it is otherwise quite slow because of the large file size
                break

        # === Create Dataframe ===
        df = pd.DataFrame(data={sensor_name: sensor_data_value, "timestamp": sensor_data_timestamp}) #convert to a dataframe
        df = df.astype({sensor_name: float, "timestamp": "datetime64[ns]"}) # set datatypes
        df.set_index("timestamp", inplace=True) # set the timestamp as the index
        df.sort_index(inplace=True) # sort the dataframe
    
        data[sensor_name] = df # add the sensor data to the data dictionary

    return data

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