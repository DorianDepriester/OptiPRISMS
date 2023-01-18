from string import Template


def fileparser(inputfname, outputfname, d):
    with open(inputfname, 'r') as f:
        src = Template(f.read())
        try:
            result = src.substitute(d)
        except KeyError as e:
            raise KeyError('Key {} in {} is undefined.'.format(e.args[0], inputfname))
    text_file = open(outputfname, "w")
    text_file.write(result)
    text_file.close()  


def CfgGenerator(d, config):
    fname = "_".join(map(str, d.values()))
    LH_name = "LHratios_{}.txt".format(fname)
    prm_name = fname + ".prm"
    d['results'] = fname
    d['LH_name'] = LH_name
    
    prisms_templates = config['PRISMS']
    fileparser(prisms_templates['prm file'], prm_name, d)
    fileparser(prisms_templates['latent hardening ratio'], LH_name, d)
    return prm_name, LH_name, fname
    
