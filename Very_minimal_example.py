from OptiPRISMS import runOptim

if __name__ == "__main__":
    # That's all. Everything you need to parametrize is actually in Config.ini
    res = runOptim(config_file='Config.ini')
    print(res)