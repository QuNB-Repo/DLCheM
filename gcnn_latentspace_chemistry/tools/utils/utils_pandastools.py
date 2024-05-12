
import pandas as pd


def filter_file_two_categories_OR(data_filepath,category_filter_1,category_filter_2,output_filepath,column_name1,column_name2):
    # Read the CSV file
    df = pd.read_csv(data_filepath)
    # Filter rows based on the value in the 'Category' column
    filtered_df = df[(df[column_name1] == category_filter_2) | (df[column_name2] == category_filter_1) ]

    # Save the filtered rows to a new CSV file
    filtered_df.to_csv(output_filepath, index=False)

def filter_file_two_categories_AND(data_filepath,category_filter_1,category_filter_2,output_filepath,column_name1,column_name2):
    # Read the CSV file
    df = pd.read_csv(data_filepath)
    # Filter rows based on the value in the 'Category' column
    filtered_df = df[(df[column_name1] == category_filter_2) & (df[column_name2] == category_filter_1) ]

    # Save the filtered rows to a new CSV file
    filtered_df.to_csv(output_filepath, index=False)


def filter_file_one_category(data_filepath,category_filter_1,output_filepath,column_name):
    # Read the CSV file
    df = pd.read_csv(data_filepath)
    filtered_df = df[(df[column_name] == category_filter_1) ]

    # Save the filtered rows to a new CSV file
    filtered_df.to_csv(output_filepath, index=False)