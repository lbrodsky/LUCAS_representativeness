def latest_version(dirs):
    for cntr in dirs:
        if not cntr.is_dir():
            continue
        versions = cntr.glob('v*')
        v_max = 0
        for v in versions:
            v_num = int(v.name[1:])
            v_max = max(v_num, v_max)

    return v_max
