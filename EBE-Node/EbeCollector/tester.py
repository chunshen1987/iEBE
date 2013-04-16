
import EbeCollector

reader = EbeCollector.EbeDBReader("testDB/collected.db")

res = reader.evaluateExpression("ecc_3(ed)")
print(res)
