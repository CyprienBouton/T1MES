from ast import literal_eval

def load_T2dash_map(value):
    try:
        # Try converting to float
        T2dash = float(value)
        return T2dash
    except ValueError:
        # Not a float? Assume it's a .npy file path
        if os.path.isfile(value) and value.endswith('.npy'):
            T2dash = np.load(value)
        else:
            raise ValueError(f"T2′ must be a .npy file or a float. Found: {value}")

def parse_key_value(arg):
    try:
        key, value = arg.split('=')
        return key, literal_eval(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Arguments must be in key=value format: got '{arg}'")

def parse_dict(pairs):
    return dict(parse_key_value(pair) for pair in pairs)