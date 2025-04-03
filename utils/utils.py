def latest_version(dirs):
    for cntr in dirs:
        versions = cntr.glob('v*')
        v_max = 0
        for v in versions:
            v_num = int(v.name[1:])
            if v_num > v_max: v_max = v_num

    return v_max
