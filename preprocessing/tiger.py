import pandas as pd
from shapely.geometry import Polygon, MultiPolygon
from shapely import wkt

# Read CSV data with specified column names
def read_data(file_path, col_names):
    return pd.read_csv(file_path, on_bad_lines='warn', names=col_names, delimiter='\t')

# Check if a given string represents a polygon
def is_polygon(str_):
    return str_.startswith("POLYGON")

def is_simple_polygon(geom):
    if isinstance(geom, Polygon):
        return not geom.interiors
    elif isinstance(geom, MultiPolygon):
        return all(not poly.interiors for poly in geom)
    return False

# Filter only the rows with valid polygon WKT
# def filter_only_simple_polygons(data):
#     return data[data.apply(is_polygon) and data.apply(is_simple_polygon)]

def filter_only_simple_polygons(data):
    # Convert WKT strings to geometries
    data['geometry'] = data['WKT'].apply(wkt.loads)
    polygons = data[data['WKT'].apply(is_polygon)]
    simple_polygons = polygons[polygons['geometry'].apply(is_simple_polygon)]

    return simple_polygons['WKT']

# Preprocess data to retain only valid polygons and save to a CSV file
def preprocess_save_polygons(in_file_path, out_file_path, col_names):
    data = read_data(in_file_path, col_names)
    valid_polygons  = filter_only_simple_polygons(data)
    save_data(valid_polygons , out_file_path)

def save_data(data, file_path):
    data.to_csv(file_path, index=False, header=False)
    
def preprocess_tiger_dataset():
    t1_input_file_path = '../../data/T1.csv'
    t1_output_csv_path = '../datasets/T1.csv'
    t1_col_names = ['WKT', 'STATEFP', 'ANSICODE', 'AREAID', 'FULLNAME', 'MTFCC', 'ALAND', 'AWATER', 'INTPTLAT', 'INTPTLON', 'PARTFLG']

    t2_input_file_path = '../../data/T2.csv'
    t2_output_csv_path = '../datasets/T2.csv'
    t2_col_names = ['WKT', 'ANSICODE', 'HYDROID', 'FULLNAME', 'MTFCC', 'ALAND', 'AWATER', 'INTPTLAT', 'INTPTLON']

    preprocess_save_polygons(t1_input_file_path, t1_output_csv_path, t1_col_names)
    preprocess_save_polygons(t2_input_file_path, t2_output_csv_path, t2_col_names)

def main():
    preprocess_tiger_dataset()

if __name__ == "__main__":
    main()