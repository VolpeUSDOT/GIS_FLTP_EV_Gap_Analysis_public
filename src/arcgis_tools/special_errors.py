class GeographicCoordinateSystemError(Exception):
    def __init__(self,message="Need a projected coordinate system to proceed."):
        self.message = message
        super().__init__(self.message)

class CoordinateSystemMismatchError(Exception):
    def __init__(self,message="Layers do not have the same coordinate system."):
        self.message = message
        super().__init__(self.message)

class FeatureClassDoesNotExistError(Exception):
    def __init__(self,message="The supplied feature class does not exist."):
        self.message = message
        super().__init__(self.message)


class WorkspaceDoesNotExistError(Exception):
    def __init__(self,message="The supplied workspace does not exist."):
        self.message = message
        super().__init__(self.message)