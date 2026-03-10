"""
Container class for pairing radar trace data with control metadata.

Defines the Trace object, which encapsulates a single-channel trace
along with its associated data array and control tuple.
"""
class Trace:
    """
    Pairs channel + data and CT metadata for a trace.

    Stores either raw or complex amplitude data along with associated channel number
    and control tuple (`ct`) for time/depth/sequence indexing.

    Parameters
    ----------
    channel : int
        Channel number associated with this trace.
    data : np.ndarray
        Trace data array; can be raw amplitudes or complex values.
    ct : tuple or list of tuples
        Control metadata.

    Attributes
    ----------
    channel, data, ct
        Same as parameters; stored without modification.
    """
    def __init__(self, channel, data, ct):
        # type: (int, np.ndarray, Union[tuple,List[tuple]]) -> None
        self.channel = channel
        self.data = data # type: np.ndarray
        self.ct = ct # type: tuple

