from OptiPRISMS import optimize

if __name__ == "__main__":
    # That's all. Everything you need to parametrize is actually in Config.ini
    res = optimize(config_file='Config.ini')
    print(res)
