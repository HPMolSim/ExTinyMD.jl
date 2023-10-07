export load_trajection

"""
load_trajection(filename::String) will load the trajection from the file named "filename".
Structure of data will be vector of vector of [id, type, x, y, z, vx, vy, vz].

"""
function load_trajection(filename::String)
    f = open(filename)
    trajection_list = []

    while !eof(f)
        readline(f)
        readline(f)

        trajection = []
        while (line = readline(f)) != ""
            data_string = split(line, ",")
            data_num = parse.(Float64, data_string)
            push!(trajection, data_num)
        end

        readline(f)

        push!(trajection_list, trajection)
    end

    close(f)

    return trajection_list
end