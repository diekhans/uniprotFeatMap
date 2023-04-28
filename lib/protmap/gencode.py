import pandas as pd

class GencodeMetadata:
    def __init__(self, gencodeMetaTsv):
        self.df = pd.read_table(gencodeMetaTsv, keep_default_na=False)

    def _doGetTrans(self, transId):
        # handle PAR_Y ids (e.g. ENST00000359512.8_PAR_Y)
        if transId.endswith("_PAR_Y"):
            transId = transId.split('_', 1)[0]
        transes = self.df[self.df.transcriptId == transId]
        if len(transes) == 0:
            raise Exception(f"can't find GENCODE metadata for '{transId}'")
        return transes.iloc[0]

    def getTrans(self, transId):
        try:
            return self._doGetTrans(transId)
        except Exception as ex:
            raise Exception(f"getTrans error for '{transId}'") from ex

    def getGeneId(self, transId):
        return self.getTrans(transId).geneId

    def getGeneName(self, transId):
        return self.getTrans(transId).geneName

    def getCodingTransIds(self):
        return list(self.df[self.df.transcriptClass == 'coding'].transcriptId)

    def getTransType(self, transId):
        return self.getTrans(transId).transcriptType
