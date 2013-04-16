
import EbeCollector

reader = EbeCollector.EbeDBReader("testDB/collected.db")

print(reader.evaluateExpression("ecc_3(e)", verbose=True))

print(reader.evaluateExpression("ecc_3(s)", verbose=True))
