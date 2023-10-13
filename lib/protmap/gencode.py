import pandas as pd
from collections import defaultdict
from pycbio.hgdata.psl import PslReader

class GencodeMetadata:
    def __init__(self, gencodeMetaTsv):
        self.df = pd.read_table(gencodeMetaTsv, keep_default_na=False)
        self.df.set_index("transcriptId", inplace=True, drop=False, verify_integrity=True)

    def _doGetTrans(self, transId):
        # handle PAR_Y ids (e.g. ENST00000359512.8_PAR_Y)
        if transId.endswith("_PAR_Y"):
            transId = transId.split('_', 1)[0]
        transes = self.df[self.df.index == transId]
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

class GencodeData:
    "all annotations"
    def __init__(self, gencodePsl):
        self.byTransId = defaultdict(list)   # older gencode PAR had same transcript ids
        for psl in PslReader(gencodePsl):
            self.byTransId[psl.qName].append(psl)
        self.byTransId.default_factory = None

    def get(self, transId, chrom):
        for psl in self.byTransId.get(transId, ()):
            if psl.tName == chrom:
                return psl
        raise Exception(f"transcript '{transId}' on chrom '{chrom}' not found in GENCODE annotations")
