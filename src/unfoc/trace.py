class Trace:
    """
    Pairs channel + data and CT metadata for a trace.
    Data is either raw amplitudes, or complex amplitudes
    """
    def __init__(self, channel, data, ct):
        # type: (int, np.ndarray, Union[tuple,List[tuple]]) -> None
        self.channel = channel
        self.data = data # type: np.ndarray
        self.ct = ct # type: tuple

