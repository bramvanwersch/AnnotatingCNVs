import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
import dash_table
import pandas as pd
import subprocess
from os import path, listdir, system
from dash.dependencies import Input, Output, State
from io import StringIO
from collections import OrderedDict
from sys import argv


#### innitial setup of global variables that do not change and are not supposed to be configurable
#this external stylesheet has a MIT liscence to should be free to use.
external_stylesheets =['https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/css/bootstrap.min.css']
ontologizer_button_active = True
# supress callback exceptions to make sure the setup function can run without layout warning interupting the building of
# the web interface.
app = dash.Dash(__name__, external_stylesheets=external_stylesheets, suppress_callback_exceptions=True)

graph_colors = ["#FFC20A", "#0C7BDC", "#1AFF1A", "#4B0092", "#994F00", "#006CD1", "#FEFE62", "#D35FB7", "#E1BE6A",
                "#40B0A6", \
                "#005AB5", "#DC3220", "#E66100", "#5D3A9B", "#1A85FF", "#D41159", "#FFC20A", "#0C7BDC", "#1AFF1A", \
                "#4B0092", "#994F00", "#006CD1", "#FEFE62", "#D35FB7", "#E1BE6A", "#40B0A6", "#005AB5", "#DC3220",
                "#E66100", \
                "#5D3A9B", "#1A85FF", "#D41159"]

vep_table_explanation = """
###### Using the table filter:
For sorting the table there is a filter box above each column. Avoid spaces in your searches. If you want to search with 
an OR operator use a "," to seperate the values for an AND operter use a "+". If you want to invert the filter and exclude 
all the words instead start your filter with a "!". Note the location column only supports the OR operator. Mixing the OR 
and AND operator is not recomended, and might give unexpected results. For the location column you can search using 
genome browser syntax eg. 1:500-5000 will return all variants on chromosome 1 inbetween 500 and 5000. Note: modifier or 
low impact does not neccesairily mean no or low impact as this is the interpretation of VEP. For more information about
the impact or consequences take a look at the link below:
"""

ontologizer_table_explanation= """
###### Running ontologizer:
For running ontologizer from this webinterface there are a few things you should know. To run ontologizer it is assumed
that the appropiate gene ontology, association and population files are present in the directory provided to this script. 
The set of genes against which is tested is the set currently displayed in the table and charts. Ontologizer needs about 
30 seconds to run so be patient. Currently the page doesnot display all results but only those that have a p-value below 
0.05. If you want to see all results take a look at thetable-(name you provided)-Parent-Child-Union-Benjamini-Hochberg.txt 
file. To view a graphical map of the go terms in their graph take a look at the png that was produced.
"""


def setup():
    """
    Function to setup some basic variables and the layout of the dash app.
    :param  df: pandas dataframe that holds the vep tab seperated values.
    :param unique_chromosomes: List containing all unique chromosomes in the pandas dataframe. This is pre calculated to
    prevent repeat calculations
    :param dfo: ontologizer dataframe holding the tab seperated data produced by ontologizer.
    """
    dir_list = listdir(result_dir)

    #try opening the file containing annotated cnvs
    try:
        vep_result_file = [result_dir + val for val in dir_list if val.startswith("added_vep")][0]
    except IndexError:
        raise FileNotFoundError("No vep file found. Make sure the file name starts with added_vep.")

    global df
    vep_version, vep_date, header_end = disect_header(vep_result_file)
    df = pd.read_csv(vep_result_file, sep="\t",skiprows=[i for i in range(header_end)],header=0, index_col=False)
    df = df.rename(index=str, columns={"#Uploaded_variation": "ID","cDNA_position":"cDNA position",\
                                       "CDS_position":"CDS position","Protein_position":"Protein position",\
                                       "Amino_acids":"Amino acids","Existing_variation":"Existing variation",
                                       "Feature_type":"Feature type"})

    global unique_chromosomes
    unique_chromosomes = [key for key, val in {name.split(":")[0]: "" for name in df["Location"]}.items()]

    #Adding white spaces to long columns, making sure that the dash table can multi-line those columns.
    df["Consequence"] = df["Consequence"].apply(add_whitespace)
    df["Extra"] = df["Extra"].apply(add_whitespace)
    df["INFO"] = df["INFO"].apply(add_whitespace)
	
    global dfo
    try:
        GO_result = [result_dir + val for val in dir_list if val.startswith("table-all_deletion_coding_cnvs")][0]
        dfo = pd.read_csv(GO_result, sep="\t", header=0)
    except IndexError:
        print("WARNING no ontologizer file found.")
        # create an empty data frame when no ontologizer file is given.
        dfo = pd.DataFrame()
    print("Loaded vep and ontologizer dataframes...")

    ##The layout of the page including the table graphs and text.
    app.layout = html.Div([
        #title
        html.H1(children='Annotation summary:'),
        html.P(children='Created from VEP {} output on {}'.format(vep_version,vep_date)),
        html.Div([
            html.Br(),
            #first text block
            dcc.Markdown(children=vep_table_explanation,
                         style={'marginRight': 20}),
            html.A('VEP table', href='https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html', target='_blank')
        ]),
        #purely to dissable the ontologizer button. Does not serve any other function
        html.Div(id='trigger', children=0, style=dict(display='none')),
        # the vep data table
        html.Div([
            html.Br(),
            dash_table.DataTable(
                id='table-sorting-filtering-graph',
                style_cell={
                'minWidth': '65px', 'maxWidth': '180px',
                'whiteSpace': 'normal',
                'textAlign': 'left'
                },
                style_cell_conditional=[
                {'if': {'column_id': 'Consequence'},
                 'width': '300px'},
                {'if': {'column_id': 'Extra'},
                 'width': '300px'},
                ],
                columns=[{'name': i, 'id': i} for i in df.columns],
                fixed_rows={'headers': True, 'data': 0},
                page_current= 0,
                page_size= 100,
                page_action='custom',

                filter_action='custom',
                filter_query='',

                sort_action='custom',
                sort_by=[],
            )
        ]),
        #components for saving csv file or gene identifiers
        html.Div([
            html.Div(dcc.Input(
                id='vep-csv-input-name',
                placeholder='Enter file name here.',
                type='text',
                ),
            style = {'margin': 'left','marginLeft': 15, 'marginTop': 20, 'column': 1}
                ),
            html.Div(html.Button(
                children='Save file',
                id='vep-create-csv-button'
                ),
                style={'margin': 'left','marginLeft': 10, 'marginTop': 20, 'column': 2}
            ),
            html.Div(dcc.RadioItems(
                id="radio-choose-download",
                options=[
                    {'label': 'Full table', 'value': 'full'},
                    {'label': 'Gene identifiers', 'value': 'ids'},
                ],
                value='full',
                labelStyle={'display': 'inline-block'}
            ),
                style={'margin': 'left', 'marginLeft': 10,'marginRight': 10, 'marginTop': 23, 'column': 3}
            ),
            html.Label(
                id='vep-output-label',
                children='',
                style={'margin': 'left', 'marginLeft': 10, 'marginTop': 20, 'column': 4}
            ),
        ],
        className = 'row'
        ),
        #containers for the graphs to go into.
        html.Div(id="type-graph-container",
                 className='row'
                 ),
        html.Div(id="consequence-graph-container",
                 className='row'
                 ),
        html.Div(id="chromosome-graph-container",
                 className='row'
                 ),
        #second text block
        html.Div([
            html.Br(),
            dcc.Markdown(children=ontologizer_table_explanation,
                         style={'marginRight': 20}),
        ]),
        #input for running the ontologizer and holding the table
        html.Div([html.Div([
                html.Div(dcc.Input(
                        id='ontologizer-run-input-name',
                        placeholder='Enter file name here.',
                        type='text',
                        ),
                        style={'margin': 'left', 'marginLeft': 15, 'marginTop': 20, 'column': 1}
                    ),
                    html.Div(html.Button(
                        id="run-ontologizer-button",
                        children="Run ontologizer",
                        ),
                        style={'margin': 'left', 'marginLeft': 10, 'marginTop': 20, 'column': 2,}
                    ),html.Label(
                        id='ontologizer-run-output-label',
                        children='',
                        style={'margin': 'left', 'marginLeft': 10, 'marginTop': 20, 'column': 3}
                    ),
                ],
                    className='row'),
                html.Div(
                    id="ontologizer-table-container"
                )]
                    ),
        #For saving the csv file of hte ontologizer output.
        html.Div([
            html.Div(dcc.Input(
                id='ontologizer-csv-input-name',
                placeholder='Enter file name here.',
                type='text',
            ),
                style={'margin': 'left', 'marginLeft': 15, 'marginTop': 20, 'column': 1}
            ),
            html.Div(html.Button(
                children='Create csv file',
                id='ontologizer-create-csv-button'
            ),
                style={'margin': 'left', 'marginLeft': 10, 'marginTop': 20, 'column': 2}
            ),
            html.Label(
                id='ontologizer-output-label',
                children='',
                style={'margin': 'left', 'marginLeft': 10, 'marginTop': 20, 'column': 3}
            ),
        ],
            className='row'
        ),
        ],
        style={'marginLeft':25}
    )
    return

def get_gene_identifiers(sort_by, filter, input_name):
    """
    Function for saving a list of gene identifiers that is currently being displayed in the table.
    :param sort_by: Sorting filter as suplied by dash.
    :param filter: Filter as suplied by dash
    :param input_name: The name specified by the user.
    :return: A string telling what happened to the file if it was saved or something else.
    """
    dataframe = filter_sort(sort_by, filter)
    #getting all unique gene identifiers
    all_gene_names = list(set([val for val in dataframe["Gene"]]))
    if input_name:
        file_loc = "{}ident_{}.txt".format(result_dir, input_name)
        if not path.exists(file_loc):
            try:
                with open(file_loc , "w") as f:
                    for elem in all_gene_names:
                        f.write(elem)
                        f.write("\n")
                return "Gene identifier file created at {}".format(file_loc)
            except OSError:
                return "Invalid file name"
        else:
            return "A file already exists with name {}".format(input_name)
    else:
        return "Please enter a name"


@app.callback(Output('run-ontologizer-button', 'disabled'),
              [Input('run-ontologizer-button', 'n_clicks'),
               Input('trigger', 'children')])
def change_button_state(nc, trigger):
    """
    Function that dissables the run-ontologizer button while ontologizer is running
    :param nc: number of from the run-ontologizer button
    :param trigger: A dummy Div component that is updated to signify when ontologizer is done running
    :return: Boolean that tells the button to be dissables (True) or enabled (False).
    """
    context = dash.callback_context.triggered[0]['prop_id'].split('.')[0]
    if context == "run-ontologizer-button":
        if nc is not None:
            return True
        else:
            return False
    else:
        return False

@app.callback(
    [Output('ontologizer-table-container', 'children'),
     Output('ontologizer-run-output-label', 'children'),
     Output('trigger', 'children')],
    [Input('run-ontologizer-button', 'n_clicks')],
    [State('table-sorting-filtering-graph', 'sort_by'),
     State('table-sorting-filtering-graph', 'filter_query'),
     State('ontologizer-run-input-name', 'value')])
def update_ontologizer(nc, sort_by, filter, input_name):
    """
    Function that runs ontologizer by executing a command line command
    :param nc: Number of clicks of the run-ontologizer-button
    :param sort_by: Sorting filter as suplied by dash.
    :param filter: Filter as suplied by dash
    :param input_name: The name specified by the user.
    :return: A list containing html and dash components that are either a label or a datatable depending on the
    ontologizer file and run.
    """
    message = ""
    if nc is not None:
        #write a file of all the gene identifiers that are currently in the table.
        message = get_gene_identifiers(sort_by, filter, input_name)
        if "file created at" in message:
            # filterign out the correct files for the analasys
            #TODO: make sure the files are named correctly so they can be found.
            dir_list = listdir(result_dir)
            global dfo
            for val in dir_list:
                if val.endswith(".obo"): obo_file = result_dir + val
                elif val.startswith("association"): association_file = result_dir + val #--> THIS FILE NEEDS TO BE NAMES association at the start
                elif val.startswith("population"): population_file = result_dir + val #--> THIS FILE NEEDS TO BE NAMED population at the start. --> both rename while runnign the pipeline.
                elif val.startswith("ident_"+input_name):study_set = result_dir + val

            try:
				#######CATCH THE FUCKING ERRORS MADE BY NEXTFLOW CRASHING -->
                command_str = "nextflow run get_go_terms.nf --obo {} --association {} --population {} --study {} --output_dir {}"\
                    .format(obo_file, association_file, population_file, study_set, result_dir)
                #still need to check if this works properly but executing using the command line is a massive pain with the localy hosted server.
                system(command_str)
                GO_file = [result_dir + val for val in listdir(result_dir) if val.startswith("table-ident_{}".format(input_name))][0]
                # recreate the message so it displays the ontologizer file created.
                message = "\\".join(GO_file.split("\\")[:-1]) + GO_file
                dfo = pd.read_csv(GO_file, sep="\t", header=0)
            except UnboundLocalError:
                message = "Cannot find the obo and/or population and association file. Please make sure they are in {}"\
                    .format(result_dir)
    if not dfo.empty:
        dffo = dfo[dfo["p.adjusted"] <= 0.05]
        if dffo.empty:
            return [[html.Br(),
                    html.Label(
                    children='No results with a p-value under the 0.05.',
                    style={'marginTop': 20}
                    )], message, nc]
        else:
            return [[html.Br(),
                    dash_table.DataTable(
                    id="ontologizer-table",
                    style_cell={
                    'minWidth': '65px', 'maxWidth': '180px',
                    'whiteSpace': 'normal',
                    'textAlign': 'left'
                    },
                    columns=[{'name': i, 'id': i} for i in dfo.columns],
                    data=dffo.to_dict('records'),
                    fixed_rows={'headers': True, 'data': 0},
                    page_current= 0,
                    page_size= 50,
                    page_action='custom',
                    )], message, nc]
    else:
        return [[html.Br(),
                html.Label(
                children='No matchig ontologizer file found.',
                style={'marginTop': 20}
                )], message, nc]

@app.callback(
    Output('vep-output-label', 'children'),
    [Input('vep-create-csv-button', 'n_clicks')],
    [State('table-sorting-filtering-graph', 'sort_by'),
     State('table-sorting-filtering-graph', 'filter_query'),
     State('vep-csv-input-name', 'value'),
     State('radio-choose-download', 'value')])
def write_from_vep_table(nc, sort_by, filter, input_name, down_choice):
    if down_choice == "full":
        return write_csv_file(nc, sort_by, filter, input_name, "vep")
    elif down_choice == "ids":
        return get_gene_identifiers(sort_by, filter, input_name)

@app.callback(
    Output('ontologizer-output-label', 'children'),
    [Input('ontologizer-create-csv-button', 'n_clicks')],
    [State('table-sorting-filtering-graph', 'sort_by'),
     State('table-sorting-filtering-graph', 'filter_query'),
     State('ontologizer-csv-input-name', 'value')])
def write_ontologizer_table(nc, sort_by, filter, input_name):
    return write_csv_file(nc, sort_by, filter, input_name, "ont")


def write_csv_file(nc, sort_by, filter, input_name, csv_name):
    #ensuring a file does not get made at starting up the program only after the button was clicked.
    if nc is not None:
        # if name is none or empty dont make te file.
        if input_name:
            if csv_name == "vep":
                dataframe = filter_sort(sort_by, filter)
            elif csv_name == "ont":
                dataframe = dfo[dfo["p.adjusted"] <= 0.05]
            file_loc = "{}\\{}.csv".format(result_dir, input_name)
            if not path.exists(file_loc):
                try:
                    dataframe.to_csv(path_or_buf=file_loc)
                    return "CSV file createt at {}".format(file_loc)
                except OSError:
                    return "Invalid file name"
            else:
                return "A file already exists with name {}".format(input_name)
        else:
            return "Please enter a name"
    return ''


##Main callback function that will call all the graphs and tables that have to be updated.
@app.callback(
    [Output('table-sorting-filtering-graph', 'data'),
     Output('type-graph-container', "children"),
     Output('consequence-graph-container', "children"),
     Output('chromosome-graph-container', "children")],
    [Input('table-sorting-filtering-graph', "page_current"),
     Input('table-sorting-filtering-graph', "page_size"),
     Input('table-sorting-filtering-graph', 'sort_by'),
     Input('table-sorting-filtering-graph', 'filter_query')])
def update_tables_graphs(page_current, page_size, sort_by, filter):
    return_list = []
    dff = filter_sort(sort_by, filter)
    page = page_current
    size = page_size
    return_list.append(dff.iloc[page * size: (page + 1) * size].to_dict('records'))
    print("Table sorted en filtered")
    return_list.append(update_graph_types(dff))
    print("Types graph generated")
    return_list.append(update_graph_consequences(dff))
    print("Consequences graph generated")
    return_list.append(update_graph_chromosome_location(dff))
    print("Chromosome graph generated")
    return return_list

def update_graph_types(dff):
    column_names = ["insertion", "deletion", "dispersed duplication", "tandem duplication"]
    # list of dictionaries containing only the unique IDs for counting of data occurances.
    type_numbers = count_types(dff)
    try:
        type_percent = [round(x/sum(type_numbers) *100, 2) for x in type_numbers]
        data = get_ordered_dict(["type", "amount", "percent"], list(zip(type_numbers, type_percent)), column_names)
    except ZeroDivisionError:
        data = OrderedDict([["types",[0]],["amount",[0]], ["percent",[0]]])
    return [html.Div(
        #pie chart displaying the amount of different types of cnvs.
            dcc.Graph(
                id='type graphs',
                className='five columns',
                figure={
                    'data': [
                        {
                            "type":"pie",
                            "labels":column_names,
                            "values":type_numbers,
                            "textinfo":'none',
                            "marker":{'colors':graph_colors,'line':{'color':'none', 'width':0}}
                        }
                    ],
                    'layout':{
                        'title':'Types of CNVs in filtered data set.',
                        }
                    }
                ),
            style={'margin': 'auto','column':1, 'width': '50%'}

    ),

        html.Div(
            dash_table.DataTable(
                data=pd.DataFrame(data).to_dict('records'),
                columns=[{'id': c, 'name': c} for c in ["type", "amount", "percent"]],
                style_cell={'textAlign': 'left'},
            ),
            style={'margin': 'auto', 'column':2, 'width': '50%'}
        )
    ]

def update_graph_consequences(dff):
    con_names, con_vals = count_consequences(dff)
    try:
        con_percent = [round(x / sum(con_vals) * 100, 2) for x in con_vals]
        data = get_ordered_dict(["consequence","amount", "percent"], list(zip(con_vals, con_percent)), con_names)
    except ZeroDivisionError:
        data = OrderedDict([["consequence",[0]],["amount",[0]], ["percent",[0]]])
    colnumcount = 2 if len(con_names)/10 > 1 else 1
    return [html.Div(
            dcc.Graph(
                id='Consequence representation',
                figure={
                    'data': [
                        {
                            "type":"pie",
                            "labels":con_names,
                            "values":con_vals,
                            "textinfo": 'none',
                            "marker":{'colors': graph_colors, 'line': {'color': 'none', 'width': 0}}
                        }
                    ],
                    'layout': {
                        'title': 'Types of consequences in filtered data set.',
                    }
                }
            ),
            style={'margin': 'auto','column':1}
    ),
        html.Div(
            dash_table.DataTable(
                data=pd.DataFrame(data).to_dict('records'),
                columns=[{'id': c, 'name': c} for c in ["consequence", "amount", "percent"]],
                style_cell={'textAlign': 'left'},
            ),
            style={'margin': 'auto','column': 2, 'columnCount':colnumcount,'width': '50%'}
        )
    ]

def update_graph_chromosome_location(dff):
    chromosomes, chrom_vals = get_chromosome_values(dff)
    try:
        chrom_percent = [round(x / sum(chrom_vals) * 100, 2) for x in chrom_vals]
        data = get_ordered_dict(["Chromosomes","Amount", "Percent"], list(zip(chrom_vals, chrom_percent)), chromosomes)
    except ZeroDivisionError:
        data = OrderedDict([["Chromosomes",[0]],["Amount",[0]], ["Percent",[0]]])
    colnumcount = 2 if len(chromosomes)/10 > 1 else 1
    return [html.Div(
        dcc.Graph(
            id="chromosome graph",
            figure={
                'data': [
                    {
                    "type":"bar",
                    "x": chromosomes,
                    "y": chrom_vals,
                    }
                ],
                'layout':{
                    'title': 'CNVs by chromosome.',
                    'xaxis':{
                        "type": 'category'
                    }
                }
            }
        ),
        style={'margin': 'auto', 'column': 1}
    ),
        html.Div(
            dash_table.DataTable(
                data=pd.DataFrame(data).to_dict('records'),
                columns=[{'id': c, 'name': c} for c in ["Chromosomes", "Amount", "Percent"]],
                style_cell={'textAlign': 'left'},
            ),
            style={'margin': 'auto', 'column': 2, 'columnCount':colnumcount, 'width': '50%'}
        )
    ]


#function that filters and sorts the data based on a given sort and filter query.
#functuion be called over and over again to sort the data. Not how fast this is and how to circumvent it.
def filter_sort(sort_by, filter):
    filtering_expressions = filter.split(' && ')
    dff = df
    for filter_part in filtering_expressions:
        col_name, filter_value = split_filter_part(filter_part)
        if col_name == "Location":
            valid_rows = []
            for loc in filter_value.split(","):
                valid_rows += filter_locations(dff, loc, col_name)
            dff = dff.loc[valid_rows]
        elif col_name is not None:
            try:
                # create regex pattern with options seperated by comma treated as or and with a plus as and. This system is not
                # perfect and mixign the two can result in unexpected results.
                # add an exclusion character ! if infront of value the value is not matched
                if filter_value.startswith("!"):
                    match = False
                    filter_value = filter_value[1:]
                else:
                    match = True
                filter_values = "|".join(filter_value.split(","))
                filter_values = "".join(["(?=.*{})".format(x) for x in filter_values.split("+")])
            except AttributeError:
                # skip if there is no filter --> for none type values.
                continue
            dff = dff.loc[dff[col_name].str.contains(filter_values) == match]
    #NOTE: multiple sorts will overwrite each other. sorting inplace does not make allot of sense because the order is
    #not apparent
    if len(sort_by):
        for col in sort_by:
            col_name = col["column_id"]
            if col_name == "Location":
                col_values = list(set([val for val in dff[col_name]]))
                if col["direction"] == "asc":
                    col_values.sort(key = lambda x: sort_by_chromosome(x))
                    sort_categories =col_values
                else:
                    col_values.sort(key = lambda x: sort_by_chromosome(x))
                    sort_categories = col_values[::-1]
                dff["Location"] = pd.Categorical(dff["Location"], sort_categories)
                dff = dff.sort_values(["Location"])
            elif col_name == "ID":
                col_values = list(set([val for val in dff[col_name]]))
                if col["direction"] == "asc":
                    col_values.sort(key = lambda x: sort_by_IDs(x))
                    sort_categories =col_values
                else:
                    col_values.sort(key = lambda x: sort_by_IDs(x))
                    sort_categories =col_values[::-1]
                dff["ID"] = pd.Categorical(dff["ID"], sort_categories)
                dff.sort_values(["ID"])
#if not location or ID sort like it normaly would.
            else:
                dff = dff.sort_values(col['column_id'],ascending=[col['direction'] == 'asc'], inplace=False)
    return dff


def sort_by_chromosome(choromosome_loc):
    """
    Function that defines a tuple for sorting chromosome locations
    :param choromosome_loc: String containing the location of the cnv in the chromosme in a genome browser syntax (1:5-10)
    :return: Tuple containing in this order a boolean determining if the chromosme is purely a number. The value of the
    chrom itself that is an integer or string depending on the previous value. Integer respresenting the start of the
    cnv and a integer representing the end.
    """
    chrom = choromosome_loc.split(":")[0]
    if "-" in choromosome_loc:
        start, end = choromosome_loc.split(":")[1].split("-")
    # in case of insertions
    else:
        start = choromosome_loc.split(":")[1]
        end = start
    # see if the chromosome can be made into an integer If this is the
    # case False is put at the start to make the numbered chromosomes come first
    try:
        return (False, int(chrom), int(start), int(end))
    except ValueError:
        return (True, chrom, int(start), int(end))

def sort_by_IDs(ID):
    """
    Function that defines a tuple for sorting IDs.
    :param ID: String representing the ID of the given cnv
    :return: Tuple containing the numerical part of the id and a boolean representing if a dot i was in the ID.
    """
    ID_parts = ID.split(".")
    if len(ID_parts) == 2:
        return (ID_parts[0], True)
    return (ID_parts[0], False)

def filter_locations(dff, filter_value, col_name):
    valid_rows = []
    try:
        chrom, startend = str(filter_value).split(":")
        start, end = startend.split("-")
        #check if chrom only contains numbers

        dff = dff.loc[dff[col_name].str.contains("^{}".format(chrom))]
        for key, col in dff[col_name].items():
            if "-" in col:
                startc, endc = col.split("-")
                startc = startc.split(":")[1]
            else:
                startc = col.split(":")[1]
                endc = startc
            if int(startc) >= int(start) and int(endc) <= int(end):
                valid_rows.append(key)
    except ValueError:
        # skip if there is no filter --> for none type values.
        return valid_rows
    return valid_rows

def split_filter_part(filter_part):
    if "contains" in filter_part:
        name_part, value_part = filter_part.split("contains", 1)
        name = name_part[name_part.find('{') + 1: name_part.rfind('}')]
        value = value_part.strip()
        v0 = value_part[0]
        if (v0 == value_part[-1] and v0 in ("'", '"', '`')):
            value = value[1: -1].replace('\\' + v0, v0)
        return name, value
    return [None] * 2


# functions for cal;culating statistics in graphs
def count_types(dff):
    #only get unique IDs
    dff = dff.drop_duplicates('ID')
    dff = dff[~dff["INFO"].str.contains("INS:DISPERSED")]
    return [len(dff[dff["Allele"] == "insertion"]),
            len(dff[dff["Allele"] == "deletion"]),
            len(dff[dff ["INFO"].str.contains("DUP:DISPERSED")]),
            len(dff[dff["INFO"].str.contains("DUP:TANDEM")])]

def count_consequences(dff):
    con_list = ['transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained', 'frameshift_variant',
     'stop_lost', 'start_lost', 'transcript_amplification', 'inframe_insertion', 'inframe_deletion', 'missense_variant',
     'protein_altering_variant', 'splice_region_variant', 'incomplete_terminal_codon_variant',
     'start_retained_variant', 'stop_retained_variant', 'synonymous_variant', 'coding_sequence_variant',
     'mature_miRNA_variant', '5_prime_UTR_variant', '3_prime_UTR_variant', 'non_coding_transcript_exon_variant',
     'intron_variant', 'NMD_transcript_variant', 'non_coding_transcript_variant', 'upstream_gene_variant',
     'downstream_gene_variant', 'TFBS_ablation', 'TFBS_amplification', 'TF_binding_site_variant',
     'regulatory_region_ablation', 'regulatory_region_amplification', 'feature_elongation', 'regulatory_region_variant',
     'feature_truncation', 'intergenic_variant']
    val_dict = {}
    for con in con_list:
        val_dict[con] = len(dff[dff["Consequence"].str.contains(con)])
    return [name for name, amnt in val_dict.items() if amnt != 0], [amnt for name, amnt in val_dict.items() if amnt != 0]

def get_chromosome_values(dff):
    val_dict = {}
    dff = dff.drop_duplicates('ID')
    dff = dff[~dff["INFO"].str.contains("INS:DISPERSED")]
    for chrom in unique_chromosomes:
        val_dict[chrom] = len(dff[dff["Location"].str.contains("^{}".format(chrom))])
    return [str(name) for name, amnt in val_dict.items() if amnt != 0], [amnt for name, amnt in val_dict.items() if amnt != 0]

def get_ordered_dict(column_names, data_rows, row_names):
    ret_list = [[name, []] for name in column_names[1:]]
    for x in range(len(data_rows)):
        for y in range(len(column_names[1:])):
            ret_list[y][1].append(data_rows[x][y])
    ret_list.insert(0, [column_names[0], row_names])
    return OrderedDict(ret_list)

def add_whitespace(column):
    column = ", ".join(column.split(","))
    column = "; ".join(column.split(";"))
    return column

def disect_header(myfile):
    header = []
    with open(myfile) as f:
        for line in f:
            if line.startswith("##"):
                header.append(line)
            else:
                break
    vep_version = header[0].split("PREDICTOR ")[1]
    vep_date = header[1].split("at ")[1]
    return vep_version, vep_date, len(header)

if __name__ == '__main__':
    global result_dir
    result_dir = argv[1]
    if not path.exists(result_dir):
        raise FileNotFoundError("Cannot locate directory: {}. Make sure the directory exists".format(result_dir))
    setup()
    app.run_server(debug=True)
