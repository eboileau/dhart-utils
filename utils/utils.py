""" Helper functions.
"""

import sys

import logging
logger = logging.getLogger(__name__)


# logging_utils


def add_logging_options(parser, default_log_file=""):
    
    """\
    This function add options for logging to an argument parser. In 
    particular, it adds options for logging to a file, stdout and stderr.
    In addition, it adds options for controlling the logging level of each
    of the loggers, and a general option for controlling all of the loggers.

    Parameters
    ----------
    parser
        an argument parser
    """

    logging_options = parser.add_argument_group("logging options")

    default_log_file = ""
    logging_level_choices = ['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']
    default_logging_level = 'WARNING'
    default_specific_logging_level = 'NOTSET'

    logging_options.add_argument('--log-file', help="This option specifies a file to "
        "which logging statements will be written (in addition to stdout and "
        "stderr, if specified)", default=default_log_file)
    logging_options.add_argument('--log-stdout', help="If this flag is present, then "
        "logging statements will be written to stdout (in addition to a file "
        "and stderr, if specified)", action='store_true')
    logging_options.add_argument('--no-log-stderr', help="Unless this flag is present, then "
        "logging statements will be written to stderr (in addition to a file "
        "and stdout, if specified)", action='store_true')

    logging_options.add_argument('--logging-level', help="If this value is specified, "
        "then it will be used for all logs", choices=logging_level_choices,
        default=default_logging_level)
    logging_options.add_argument('--file-logging-level', help="The logging level to be "
        "used for the log file, if specified. This option overrides "
        "--logging-level.", choices=logging_level_choices, 
        default=default_specific_logging_level)
    logging_options.add_argument('--stdout-logging-level', help="The logging level to be "
        "used for the stdout log, if specified. This option overrides "
        "--logging-level.", choices=logging_level_choices, 
        default=default_specific_logging_level)
    logging_options.add_argument('--stderr-logging-level', help="The logging level to be "
        "used for the stderr log, if specified. This option overrides "
        "--logging-level.", choices=logging_level_choices, 
        default=default_specific_logging_level)


def update_logging(args, logger=None, 
        format_str='%(levelname)-8s %(name)-8s %(asctime)s : %(message)s'):

    """\
    This function interprets the logging options in args. Presumably, these
    were added to an argument parser using add_logging_options.

    Parameters
    ----------
    args
        a namespace with the arguments added by add_logging_options

    logger
        a logger which will be updated. If None is given, then the default
        logger will be updated.

    format_str
        The logging format string. Please see the python logging documentation
        for examples and more description.
    """

    # find the root logger if another logger is not specified
    if logger is None:
        logger = logging.getLogger('')
            
    logger.handlers = []

    # set the base logging level
    level = logging.getLevelName(args.logging_level)
    logger.setLevel(level)

    # now, check the specific loggers

    if len(args.log_file) > 0:
        h = logging.FileHandler(args.log_file)
        formatter = logging.Formatter(format_str)
        h.setFormatter(formatter)
        if args.file_logging_level != 'NOTSET':
            l = logging.getLevelName(args.file_logging_level)
            h.setLevel(l)
        logger.addHandler(h)

    if args.log_stdout:
        h = logging.StreamHandler(sys.stdout)
        formatter = logging.Formatter(format_str)
        h.setFormatter(formatter)
        if args.stdout_logging_level != 'NOTSET':
            l = logging.getLevelName(args.stdout_logging_level)
            h.setLevel(l)
        logger.addHandler(h)

    log_stderr = not args.no_log_stderr
    if log_stderr:
        h = logging.StreamHandler(sys.stderr)
        formatter = logging.Formatter(format_str)
        h.setFormatter(formatter)
        if args.stderr_logging_level != 'NOTSET':
            l = logging.getLevelName(args.stderr_logging_level)
            h.setLevel(l)
        logger.addHandler(h)
        
        
# utils    






