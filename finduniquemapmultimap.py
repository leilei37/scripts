#file_name_template = '%dLog.final.out'
lines_out = []
with open('alllogout.txt', 'r') as f_in:
    print(f_in)
    lines1 = f_in.readlines()
    print(lines1)
    with open('uniquely_mapped_output.txt', 'w') as f_out:
        for line in lines1:
            linee = line[:-1]
            with open(linee, 'r') as f:
                lines = f.readlines()
                nine_line = lines[8]
                split_line = nine_line.split('|')
                number = split_line[1].strip()
                twentyfour_line = lines[23]
                split_line2 = twentyfour_line.split('|')
                multinumber = split_line2[1].strip()
                lines_out.append(f'{linee}\t{number}\t{multinumber}')
        for line_out in lines_out:
            f_out.write(line_out+'\n')
        f_out.write('\n')