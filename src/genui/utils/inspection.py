import sys
import re
import pkgutil
import importlib
import inspect
import qsprpred
import sklearn
from sklearn.base import RegressorMixin, ClassifierMixin
from sklearn.utils._param_validation import Interval, StrOptions
from sklearn.metrics._dist_metrics import DistanceMetric
import numpy as np
from django.urls import path, include
from django.db import models
from genui import apps


def importFromPackage(package, module_name, exception=True):
    package_name = package
    if type(package) != str:
        package_name = package.__name__
    try:
        return importlib.import_module(f'{package_name}.{module_name}')
    except ModuleNotFoundError as exp:
        if exception:
            raise exp
        else:
            return None


def importModuleWithException(module, *args, message=True, throw=False, **kwargs):
    try:
        return importlib.import_module(module, *args, **kwargs)
    except ModuleNotFoundError as exp:
        if message:
            print(f"WARNING: Failed to find core modelling module: {module}. It will be skipped.")
        if throw:
            raise exp
        return


def getSubclassesFromModule(base_cls, module):
    """
    Fetch all subclasses of a given base class from a module.

    :param base_cls: The base class.
    :param module: The module.
    :return:
    """

    ret = []
    for item in inspect.getmembers(module, inspect.isclass):
        if issubclass(item[1], base_cls):
            ret.append(item[1])
    return ret


def getSubclasses(cls):
    """
    Fetch all existing subclasses of a class.

    :param cls:
    :return:
    """

    return set(cls.__subclasses__()).union([s for c in cls.__subclasses__() for s in getSubclasses(c)])


def findSubclassByID(base, module, id_attr: str, id_attr_val: str):
    """
    Function to fetch a given class from a certain module.
    It is identified by both its base class and a value of a
    specific identifying attribute of the class.

    :param base: The base class of the class we are looking for.
    :param module: Module containing the class of interest.
    :param id_attr: The name of the identifying attribute on the class.
    :param id_attr_val: The value of the searched attribute.
    :return:
    """

    all_subclasses = getSubclasses(base)

    for class_ in all_subclasses:
        name = class_.__name__
        if not hasattr(module, name):
            continue

        if not hasattr(class_, id_attr):
            raise Exception(
                "Unspecified ID attribute on a class where it should be defined. Check if the class is properly annotated: ",
                repr(class_))

        if not (id_attr_val == getattr(class_, id_attr)):
            continue

        return class_

    raise LookupError(f'Could not find any valid subclass of {base.__name__} in module {module.__name__}.')


def getFullName(obj, moduleOnly=False):
    module = inspect.getmodule(obj)
    if module is None or module == str.__class__.__module__:
        return obj.__name__ if not moduleOnly else None
    else:
        return module.__name__ + '.' + obj.__name__ if not moduleOnly else module.__name__


def getObjectAndModuleFromFullName(name):
    module = importlib.import_module('.'.join(name.split(".")[0:-1]))
    obj = getattr(module, name.split(".")[-1])
    return obj, module


def discover_extensions(packages):
    ret = []
    for package_name in packages:
        try:
            package = importlib.import_module(package_name)
        except ModuleNotFoundError as exp:
            sys.stderr.write(f"WARNING: Failed to find extensions module: {package_name}. Skipping...\n")
            continue

        subs = package.__all__
        for sub in subs:
            ret.append(f"{package_name}.{sub}")
    return ret


def disover_app_urls_module(app_name, parent=None):
    app = importlib.import_module(app_name)

    error_skipped = False
    if not importFromPackage(app, 'genuisetup', exception=False):
        sys.stderr.write(f"WARNING: Expected to find the {app_name}.genuisetup module, but it was not found.\n")
        error_skipped = True
    elif not importFromPackage(app, 'urls', exception=False):
        # sys.stderr.write(f"WARNING: No {app_name}.urls module found for app:  {app_name}\n")
        error_skipped = True
    elif not hasattr(app.genuisetup, 'PARENT'):
        if parent == 'genui':
            return app.urls
        else:
            return None
    elif app.genuisetup.PARENT != parent:
        return None

    if error_skipped:
        # sys.stderr.write(f"WARNING: Skipping urls discovery for: {app_name}\n")
        return None

    return app.urls


def discover_apps_urls(app_names, prefix='', app_names_as_root=False):
    urls = []
    for app in app_names:
        urls_module = disover_app_urls_module(app, parent='genui')
        if urls_module:
            urls.append(path(
                f'{prefix.rstrip("/") + "/" if prefix else ""}{app.split(".")[1] + "/" if app_names_as_root else ""}',
                include(urls_module.__name__)
            )
            )

    return urls


def discover_extensions_urlpatterns(parent):
    ret = []
    for extension_name in apps.extensions():
        urls = disover_app_urls_module(extension_name, parent)
        if urls:
            ret.extend(urls.urlpatterns)
    return ret


def get_non_abstract_classes_from_module(module):
    if isinstance(module, str):
        module = importlib.import_module(module)

    classes = []
    for name, obj in inspect.getmembers(module):
        if inspect.isclass(obj) and not inspect.isabstract(obj) and obj.__module__ == module.__name__:
            classes.append(name)
    return classes


def get_sklearn_models():
    regressors = {}
    classifiers = {}
    search_modules = {"ensemble", "gaussian_process", "kernel_ridge", "linear_model", "neighbors", "neural_network",
                      "svm", "tree"}

    for _, module_name, _ in pkgutil.walk_packages(sklearn.__path__, prefix="sklearn."):
        try:
            if set(module_name.split(".")) & search_modules:
                module = importlib.import_module(module_name)

                for name, obj in inspect.getmembers(module):
                    if inspect.isclass(obj) and obj.__module__ == module_name:
                        if not inspect.isabstract(obj) and not obj.__name__.startswith("_") and not "CV" in name:
                            if issubclass(obj, RegressorMixin):
                                regressors[name] = f"{module_name}.{name}"
                            elif issubclass(obj, ClassifierMixin):
                                classifiers[name] = f"{module_name}.{name}"
        except Exception as e:
            pass

    return regressors, classifiers


def get_metrics():
    curve_metrics = ["precision_recall_curve", "roc_curve", "det_curve"]

    classification_metrics = [
        'accuracy', 'average_precision', 'balanced_accuracy', 'f1', 'jaccard','matthews_corrcoef', 'precision', 'recall',
        'roc_auc',
    ]

    regression_metrics = [
        'd2_absolute_error_score', 'explained_variance', 'neg_max_error',
        'neg_mean_absolute_error', 'neg_mean_absolute_percentage_error',
        'neg_mean_gamma_deviance', 'neg_mean_poisson_deviance', 'neg_mean_squared_error',
        'neg_mean_squared_log_error', 'neg_median_absolute_error',
        'neg_root_mean_squared_error', 'neg_root_mean_squared_log_error', 'r2'
    ]

    cls_metrics = curve_metrics + classification_metrics

    return {"classification": cls_metrics, "regression": regression_metrics}


def get_data_splits():
    qsprpred_split_classes = get_non_abstract_classes_from_module('qsprpred.data.sampling.splits')
    qsprpred_split_classes = [f'qsprpred.data.sampling.splits.{name}' for name in qsprpred_split_classes]

    sklearn_split_classes = get_non_abstract_classes_from_module('sklearn.model_selection._split')
    sklearn_split_classes = [f'sklearn.model_selection.{name}' for name in sklearn_split_classes if
                             any(suffix in name for suffix in ['Split', 'Fold', 'CV']) and not name.startswith('_')]

    return qsprpred_split_classes + sklearn_split_classes


def parse_interval(constraint):
    constraint = str(constraint)
    type_ = "int" if "int" in constraint else "float" if "float" in constraint else None
    if type_ is None:
        return None
    range_str = " ".join(constraint.rsplit(" ", 2)[1:])
    lb = range_str[0]
    rb = range_str[-1]
    range_str = range_str[1:-1]
    l, r = range_str.split(",")
    l = l.strip()
    r = r.strip()
    interval = {}
    if l != "-inf": interval["min"] = float(l)
    if r != "inf": interval["max"] = float(r)
    if "min" in interval: interval["leq"] = "bq" if lb == "[" else "b"
    if "max" in interval: interval["req"] = "sq" if rb == "]" else "s"
    value = {"type": type_, "interval": interval}
    return value


def _constraint_processor(constraint):
    if type(constraint) == Interval:
        return parse_interval(constraint)
    elif type(constraint) == StrOptions:
        return {"type": "str", "choices": [o for o in constraint.options]}
    elif inspect.isclass(constraint) and issubclass(constraint, DistanceMetric):
        return {"type": "str", "choices": ["euclidean", "manhattan", "chebyshev'"]}
    elif constraint == "boolean" or (
            inspect.isclass(constraint) and (issubclass(constraint, bool) or issubclass(constraint, np.bool_))):
        return {"type": "bool"}
    elif inspect.isclass(constraint) and (issubclass(constraint, list) or issubclass(constraint, dict)):
        return None
    else:
        return None


DISCARDED_PARAMS = ["self", "random_state", "verbose", "n_jobs", ]


def get_sklearn_params_with_constraints(module_name, class_=None):
    if class_ is None:
        class_ = module_name.split('.')[-1]
        module_name = '.'.join(module_name.split('.')[:-1])

    module = importlib.import_module(module_name)
    class_ = getattr(module, class_)
    signature = inspect.signature(class_)

    constraints = getattr(class_, "_parameter_constraints", {})

    params = {}
    for param_name, param in signature.parameters.items():
        if param_name in DISCARDED_PARAMS:
            continue
        default = None if param.default is inspect.Parameter.empty else param.default
        constraint = constraints.get(param_name, None)

        if constraint is not None:
            constraint = [_constraint_processor(c) for c in constraint if c is not None]
            constraint = [c for c in constraint if c is not None and c != False]
            if isinstance(default, str):
                constraint = [c for c in constraint if c["type"] == "str"]
            constraint = constraint[0] if len(constraint) > 0 else None

        if default is None and constraint is not None:
            c = constraint
            ctype = c["type"]
            if ctype == "int" or ctype == "float":
                default = c["interval"]["min"] if c["interval"]["leq"] != "b" else c["interval"]["min"] + 0.001
            elif ctype == "str":
                default = c["choices"][0]
            elif ctype == "bool":
                default = False
        if constraint is not None:
            params[param_name] = {
                "value": default,
                "constraint": constraint,
            }

    return params


sklearn_regressors, sklearn_classifiers = get_sklearn_models()
SKLEARN_MODELS = sklearn_regressors | sklearn_classifiers
SKLEARN_MODELS_PARAMS = {k: get_sklearn_params_with_constraints(v) for k, v in SKLEARN_MODELS.items()}
METRICS = get_metrics()
DATA_SPLITS = get_data_splits()
SCAFFOLDS = get_non_abstract_classes_from_module("qsprpred.data.chem.scaffolds")
CLUSTERING = get_non_abstract_classes_from_module("qsprpred.data.chem.clustering")


def get_default_params(class_=None, module_name=None):
    if class_ is None:
        class_ = module_name.split('.')[-1]
        module_name = '.'.join(module_name.split('.')[:-1])
    params = {}
    module = importlib.import_module(module_name)
    class_ = getattr(module, class_)
    signature = inspect.signature(class_)
    for param_name, param in signature.parameters.items():
        if param.default is inspect.Parameter.empty:
            params[param_name] = None
        else:
            if not isinstance(param.default, (int, float, str)):
                params[param_name] = str(param.default)
            else:
                params[param_name] = param.default
    return params


DISCARD_NAMES = ["polymorphic_ctype", "trainingStrategy"]


def get_default_params_django(class_=None, module_name=None):
    django_field2python = {
        "CharField": "str",
        "TextField": "str",
        "IntegerField": "int",
        "FloatField": "float",
        "BooleanField": "bool",
        "DateField": "datetime.date",
        "DateTimeField": "datetime.datetime",
    }
    if class_ is None:
        class_ = module_name.split('.')[-1]
        module_name = '.'.join(module_name.split('.')[:-1])
    module = importlib.import_module(module_name)
    model_class = getattr(module, class_)
    parameters = {}
    if not hasattr(model_class, "arguments"):
        for field in model_class._meta.get_fields():
            if isinstance(field, models.Field) and not field.auto_created:
                field_type = field.get_internal_type()
                field_type = django_field2python.get(field_type, field_type)
                default_value = field.default if field.default != models.NOT_PROVIDED else None
                if field.name not in DISCARD_NAMES:
                    parameters[field.name] = {"type": field_type, "value": default_value}
    else:
        parameters = model_class.arguments
    return parameters


def camel_to_snake(name):
    return re.sub(r'(?<!^)(?=[A-Z])', '_', name).lower()


def snake_to_camel(name):
    words = name.split('_')
    return words[0] + ''.join(word.capitalize() for word in words[1:])
