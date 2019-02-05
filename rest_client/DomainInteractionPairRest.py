import json
from rest_client.GetRest import GetRest
from rest_client.PostRest import PostRest

class DomainInteractionPairAPI(object):
    """
    This class manage the requests for the Domain Interaction Pair objects into the restAPI

    :param function: the name of the function to access in the rest API
    :type function: string
    """

    def __init__(self, function='domainintpair/'):
        """
        Initialization of the class

        :param function: name of the function

        :type function: string (url)

        """
        self.function = function


    def get_all(self):
        """
        get all the domains interaction pairs on the database

        :return: json file with all the data
        :rtype: string (json format)
        """
        result_get = GetRest(function = self.function).performRequest()
        return result_get

    def setDomainInteractionPair(self, jsonData):
        """
        set new domain interaction pair in the database

        :return: json file with the last domain created
        :rtype: string (json format)
        """
        jsonData = json.dumps(jsonData)
        result_post = PostRest(function = self.function, dataDict = jsonData).performRequest()
        return result_post