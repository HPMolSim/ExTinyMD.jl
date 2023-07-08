export read_trajection

function read_trajection(filename::String, n_atoms::T) where{T<:Integer}
    # 打开文件
    f = open(filename)

    # 初始化矩阵列表
    matrix_list = []

    # 循环读取文件
    while !eof(f)
        # 读取第一行数字
        readline(f)

        num = n_atoms
        # 读取矩阵
        matrix = zeros(num, 3)
        for i in 1:num
            line = readline(f)
            data_string = split(line, ",")
            data_num = parse.(Float64, data_string)
            matrix[i, :] = data_num
        end

        # 跳过空行
        readline(f)
        readline(f)

        # 将矩阵添加到列表中
        push!(matrix_list, matrix)
    end

    # 关闭文件
    close(f)

    # 返回矩阵列表
    return matrix_list
end
