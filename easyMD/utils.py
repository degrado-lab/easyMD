import importlib.util
import sys

def import_function_from_path(path, function_name):
    spec = importlib.util.spec_from_file_location("module.name", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["module.name"] = module
    spec.loader.exec_module(module)
    return getattr(module, function_name)

def check_required_params(params, required_params):
    for param in required_params:
        if param not in params:
            raise ValueError(f'Required parameter \'{param}\' not found in input parameters. Please provide this parameter.')