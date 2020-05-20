# You should read carefully the correct way to acknowledge the use of this data set / the
# appropriate papers to cite at the top of this readme:
# https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt

using FTPClient, CSV, DataFrames, Dates



#
# Functionality pertaining to downloading and decompressing the GHCN data.
#

"""
    const ftp_url

Location in the FTP server to connect to.
"""
const ftp_url = "ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/"

"""
    const station_metadata_fname

Name of the file containing station metadata.
"""
const station_metadata_fname = "ghcnd-stations.txt"

"""
    const countries_fname

Name of the file containing country-code -- code pairing data set. See
`download_countries` and `parse_countries_metadata` for functionality.
"""
const countries_fname = "ghcnd-countries.txt"

"""
    const inventory_fname

Name of the inventory file in the GHCN daily data.
"""
const inventory_fname = "ghcnd-inventory.txt"

"""
    const ghcnd_data

Name of the file containing all of the ghcn daily data.
"""
const ghcnd_data = "ghcnd_all"

"""
    open_ftp(url)

Open a connection with NOAA's FTP server, at the Global Historical Climate Network (GHCN)
daily data.
"""
open_ftp() = FTP(ftp_url)

# Helper function.
download_file(ftp, fname) = download(ftp, fname, fname)

"""
    download_station_metadata()

Downloads the metadata associated with all stations in the GHCN data set.
"""
function download_station_metadata()
    ftp = open_ftp()
    download_file(ftp, station_metadata_fname)
    close(ftp)
    return nothing
end

"""
    download_countries()

Downloads data that pairs specifies the two character code for each country present within
the GHCN data set.
"""
function download_countries()
    ftp = open_ftp()
    download_file(ftp, countries_fname)
    close(ftp)
    return nothing
end

"""
    download_inventory()

Downloads data explaining the variables available at each station, and the start / end date
of the available time.
"""
function download_inventory()
    ftp = open_ftp()
    download_file(ftp, inventory_fname)
    close(ftp)
    return nothing
end

"""
    download_all_data()

Downloads all of the data from the GHCN network.
"""
function download_all_data()
    ftp = open_ftp()
    download_file(ftp, ghcnd_data * ".tar.gz")
    close(ftp)
    return nothing
end

"""
    uncompress_all_data(fname="ghcnd_all.tar.gz")

`download_all_data` obtained a `tar.gz`. This function uncompresses it. Assumes a working
installation of tar.
"""
function uncompress_all_data(fname="ghcnd_all.tar.gz")
    run(`tar xzvf $(ghcnd_data * "tar.gz")`)
end



#
# Functionality pertaining to loading the data in a useful format.
#

"""
    load_station_metadata()

Loads the station metadata detailed in section 4 of [1], returning it in a `DataFrame` with
columns named according the the `Variable` column below:

```
------------------------------
Variable   Columns   Type
------------------------------
ID            1-11   Character
LATITUDE     13-20   Real
LONGITUDE    22-30   Real
ELEVATION    32-37   Real
STATE        39-40   Character
NAME         42-71   Character
GSN FLAG     73-75   Character
HCN/CRN FLAG 77-79   Character
WMO ID       81-85   Character
------------------------------
```

[1] - https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt
"""
function load_station_metadata()
    raw_data = CSV.read(station_metadata_fname; header=0, delim='\t').Column1
    N_stations = length(raw_data)

    parsed_data = (
        ID = Vector{String}(undef, N_stations),
        LATITUDE = Vector{Float32}(undef, N_stations),
        LONGITUDE = Vector{Float32}(undef, N_stations),
        ELEVATION = Vector{Float32}(undef, N_stations),
        STATE = Vector{String}(undef, N_stations),
        NAME = Vector{String}(undef, N_stations),
        GSN_FLAG = Vector{String}(undef, N_stations),
        HCN_CRN_FLAG = Vector{String}(undef, N_stations),
        WMO_ID = Vector{String}(undef, N_stations),
    )

    for (n, row) in enumerate(raw_data)
        parsed_data.ID[n] = row[1:11]
        parsed_data.LATITUDE[n] = parse(Float32, row[13:20])
        parsed_data.LONGITUDE[n] = parse(Float32, row[22:30])
        parsed_data.ELEVATION[n] = parse(Float32, row[32:37])
        parsed_data.STATE[n] = strip(row[39:40])
        parsed_data.NAME[n] = strip(row[42:71])
        parsed_data.GSN_FLAG[n] = strip(row[73:75])
        parsed_data.HCN_CRN_FLAG[n] = strip(row[77:79])
        parsed_data.WMO_ID[n] = strip(row[81:85])
    end

    return DataFrame(parsed_data)
end

"""
    load_countries_metadata()

Loads the countries metadata detailed in section 5 of [1], returning the result as a
`DataFrame` with columns `CODE` and `NAME`.

[1] - https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt
"""
function load_countries_metadata()
    raw_data = CSV.read(countries_fname; header=0, delim='\t').Column1
    N_countries = length(raw_data)

    parsed_data = (
        CODE = Vector{String}(undef, N_countries),
        NAME = Vector{String}(undef, N_countries),
    )

    for (n, row) in enumerate(raw_data)
        parsed_data.CODE[n] = row[1:2]
        parsed_data.NAME[n] = strip(row[4:end])
    end
    return DataFrame(parsed_data)
end

"""
    load_inventories()

Loads the inventory file detailed in section 7 of [1], and returns a `DataFrame` with
columns named according to the variables listed in the file-format information below.

------------------------------
Variable   Columns   Type
------------------------------
ID            1-11   Character
LATITUDE     13-20   Real
LONGITUDE    22-30   Real
ELEMENT      32-35   Character
FIRSTYEAR    37-40   Integer
LASTYEAR     42-45   Integer
------------------------------

[1] - https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt
"""
function load_inventories()
    raw_data = CSV.read(inventory_fname; header=0, delim='\t').Column1
    N_stations = length(raw_data)

    parsed_data = (
        ID = Vector{String}(undef, N_stations),
        LATITUDE = Vector{Float32}(undef, N_stations),
        LONGITUDE = Vector{Float32}(undef, N_stations),
        ELEMENT = Vector{String}(undef, N_stations),
        FIRSTYEAR = Vector{Int}(undef, N_stations),
        LASTYEAR = Vector{Int}(undef, N_stations),
    )

    for (n, row) in enumerate(raw_data)
        parsed_data.ID[n] = row[1:11]
        parsed_data.LATITUDE[n] = parse(Float32, row[13:20])
        parsed_data.LONGITUDE[n] = parse(Float32, row[22:30])
        parsed_data.ELEMENT[n] = row[32:35]
        parsed_data.FIRSTYEAR[n] = parse(Int, row[37:40])
        parsed_data.LASTYEAR[n] = parse(Int, row[42:45])
    end

    return DataFrame(parsed_data)
end



"""
    load_data_file(station_id::String)

See section 3 of [1]...

```
------------------------------
Variable   Columns   Type
------------------------------
ID            1-11   Character
YEAR         12-15   Integer
MONTH        16-17   Integer
ELEMENT      18-21   Character
VALUE1       22-26   Integer
MFLAG1       27-27   Character
QFLAG1       28-28   Character
SFLAG1       29-29   Character
VALUE2       30-34   Integer
MFLAG2       35-35   Character
QFLAG2       36-36   Character
SFLAG2       37-37   Character
  .           .          .
  .           .          .
  .           .          .
VALUE31    262-266   Integer
MFLAG31    267-267   Character
QFLAG31    268-268   Character
SFLAG31    269-269   Character
------------------------------
```

[1] - https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt
"""
function load_data_file(station_id::String)
    file_name = joinpath(ghcnd_data, station_id * ".dly")
    raw_data = CSV.read(file_name; header=0, delim='\t').Column1
    N_months = length(raw_data)

    # Construct a `NamedTuple` to store the parsed data.
    column_names = vcat(
        [:ID, :YEAR, :MONTH, :ELEMENT],
        map(1:31) do n
            [Symbol("VALUE$n"), Symbol("MFLAG$n"), Symbol("QFLAG$n"), Symbol("SFLAG$n")]
        end...,
    )
    columns = vcat(
        [
            Vector{String}(undef, N_months),
            Vector{Int}(undef, N_months),
            Vector{Int}(undef, N_months),
            Vector{String}(undef, N_months),
        ],
        map(1:31) do n
            [
                Vector{Union{Int, Missing}}(undef, N_months),
                Vector{String}(undef, N_months),
                Vector{String}(undef, N_months),
                Vector{String}(undef, N_months),
            ]
        end...,
    )

    # Specify the beginning and end location of each column of the data.
    column_locations = vcat(
        [(1, 11), (12, 15), (16, 17), (18, 21)],
        map(1:31) do n
            start_pos = 22 + 8 * (n - 1)
            [
                (start_pos, start_pos + 4),
                (start_pos + 5, start_pos + 5),
                (start_pos + 6, start_pos + 6),
                (start_pos + 7, start_pos + 7),
            ]
        end...,
    )

    # Iterate over the data.
    for (n, row) in enumerate(raw_data)
        for col in eachindex(columns)
            start_pos = first(column_locations[col])
            end_pos = last(column_locations[col])
            columns[col][n] = _to(eltype(columns[col]), row[start_pos:end_pos])
        end
    end

    return DataFrame(columns, column_names)
end

function _to(::Type{Union{Missing, T}}, s::String) where {T<:Number}
    if s == "-9999"
        return missing
    else
        return parse(T, s)
    end
end
_to(::Type{String}, s::String) = s



#
# Common operations that reformat the data.
#


"""
    convert_to_time_series(df::DataFrame)

Converts a `DataFrame` loaded via `load_data_file` from the original format, in which each
row contains a month's worth of data for a particular ELEMENT (variable type e.g. precip,
temperature, etc), to a format in which each row contains a single day's worth of data with
for a single ELEMENT type. Only the element types already present in `df` will be included
in the new `DataFrame`.
"""
function convert_to_time_series(df::DataFrame)

    # A dataframe in which each lines containing a single day _and_ a single element.
    tmp = stack(a_file, 5:size(a_file, 2))

    # Compute day-of-month.
    tmp.DAY = map(x -> parse(Int, string(x)[6:end]), tmp.variable)

    # Filter to remove invalid dates.
    filter!(row -> Dates.validargs(Date, row.YEAR, row.MONTH, row.DAY) === nothing, tmp)

    # Merge date information into a single column.
    tmp.DATE = Date.(tmp.YEAR, tmp.MONTH, tmp.DAY)

    # Merge element and variable-type into a single underscore-separated column.
    tmp.col_index = map(zip(tmp.ELEMENT, tmp.variable)) do (element, variable)
        element * "_" * string(variable)[1:5]
    end

    # Drop redundant columns.
    tmp2 = tmp[!, [:DATE, :col_index, :value]]

    # Widen the currently very tall data frame.
    tmp3 = unstack(tmp2, [:DATE], :col_index, :value)

    # Reduce type unions over `Any` columns.
    for c in 2:size(tmp3, 2)
        name = names(tmp3)[c]
        if name[6:end] == "VALUE"
            tmp3[!, name] = convert(Vector{Union{Int, Missing}}, tmp3[!, name])
        else
            tmp3[!, name] = convert(Vector{Union{String, Missing}}, tmp3[!, name])
        end
    end

    # Re-order the data in a more natural order.
    sort!(tmp3, :DATE)

    return tmp3
end


"""
    dataset_from_region(time_interval, lat_interval, lon_interval, element)

Load all of the data from a particular region and produce a 
"""
function dataset_from_region(
    time_interval::Tuple,
    lat_interval::Tuple{Real, Real},
    lon_interval::Tuple{Real, Real},
    element::String,
)
    # Load all of the inventories.
    inventories = load_inventories()

    # Filter them to exclude all stations that don't have data that lives inside the
    # specified time period.
    filter(inventories) do row
        row.ELEMENT == element &&
            lat_interval[1] < row.LATITUDE < lat_interval[2] &&
            lon_interval[1] < row.LONGITUDE < lon_interval[2] &&
            (
                Date(row.FIRSTYEAR, 01, 01) <= time_interval[1] <= Date(row.LASTYEAR, 12, 31) ||
                Date(row.FIRSTYEAR, 01, 01) <= time_interval[2] <= Date(row.LASTYEAR, 12, 31)
            )
    end

    # What format do I want to return the data in? A matrix with NaNs / Missings?
end

# Load the metadata
countries_data = CSV.read(countries_fname)
