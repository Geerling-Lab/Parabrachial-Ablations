def MouseLookup(mouse_number, columns, df=None):
    """
    Opens up Lesion Mice Overview.txt, and finds information about mouse_number
    :param mouse_number: string with mouse number you want to look up
    :param columns: list of strings with column titles you want to look up
    :param dfs: optional pandas dataframe. This function will generate its own if one doesn't exist
    :return:
    """
    mouse_number = int(mouse_number)
    import pandas as pd
    if not df:
        df = pd.read_excel(io=r"R:/Fillan/Parabrachial Ablations/Lesion Mice Overview.xlsx", sheet_name="Overview")
    df.set_index("Mouse", inplace=True)
    return list(df.loc[mouse_number, columns])
