"""For interrogating the queue and list of completed tasks"""

                                     
__author__ = "Graham Poulter"                                        
__license__ = "GPL"

from mscanner.configuration import rc
from mscanner import scorefile


def make_queue_name(timestamp, dataset):
    """Given a time.time() stamp and task name, construct queue file name"""
    import time
    ftime = time.strftime("%Y%m%d-%H%M%S", time.gmtime(timestamp))
    return rc.queue_path / ("%s-%s.txt" % (ftime, dataset))


def parse_queue_name(fname):
    """Given path to queue file, recover task name"""
    fmatch = re.match(r"\d{8}-\d{6}-(.+)\.txt", filename.basename())
    return fmatch.group(1) if fmatch is not None else None


def get_task_status(dataset, current=None):
    """Get the status of the specified task 

    @return: (status, descriptor) where status in ["current", "done", "queue",
    "notfound"]."""
    if current is not None and dataset == current.dataset:
        # Task is the one being processed
        return "current", current
    if (rc.web_report_dir / dataset / rc.report_index).exists():
        # Task is is completed, descriptor is in the output directory
        return "done", scorefile.readDescriptor(
            rc.web_report_dir / dataset / rc.report_descriptor)
    for queue_dataset, fname in zip(*list_queue()):
        if dataset == queue_dataset:
            # Task is in the queue
            return "queue", scorefile.readDescriptor(rc.queue_path / fname)
    # Task is nowhere to be found
    return "notfound", None


def get_current_task():
    """
    @return: None if the queue is empty, otherwise a Storage
    descriptor for the currently running task.

    @note: Descriptor has additional queue_length variable for total number of
    tasks in the queue."""
    listing = rc.queue_path.files("*.txt")
    if len(listing) == 0: return None
    current = scorefile.readDescriptor(listing[0])
    current.queue_length = len(listing)
    return current


def list_queue(descriptors=False):
    """
    @return: (tasks, files) with tasks in the queue and corresponding
    filenames.

    @param descriptors: If True, tasks is a list of descriptor objects. If
    False, tasks is a list of task names."""
    tasks = []
    filenames = rc.queue_path.files("*.txt")
    for filename in filenames:
        dataset = parse_queue_name(filename)
        if dataset is not None:
            if descriptors == True:
                tasks.append(scorefile.readDescriptor(filename))
            else:
                tasks.append(dataset)
    return tasks, filenames


def list_done(descriptors=False):
    """
    @return: List of completed tasks.

    @param descriptors: If True, tasks are full descriptors. If False, results
    are just the names of the tasks.
    """
    output = []
    for fpath in rc.web_report_dir.dirs():
        if (fpath / rc.report_index).exists():
            if descriptors == True:
                output.append(
                    scorefile.readDescriptor(fpath / rc.report_descriptor))
            else:
                output.append(fpath.basename())
    return output


def list_visible():
    """Tasks which should be visible in the output list
    
    @return: List of task descriptors"""
    donetasks = sorted(list_done(descriptors=True), 
                       key=lambda x:x.timestamp, reverse=True)
    return [ d for d in donetasks if "hidden" not in d or d.hidden == False ]
