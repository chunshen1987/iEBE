
import EbeCollector

reader = EbeCollector.EbeDBReader("testDB/collected.db")

reader.evaluateExpression("dN/dydpT(0.5)(pion)", verbose=True)

print(reader.evaluateExpression("ecc_3(e)", verbose=True))

print(reader.evaluateExpression("ecc_3(s)", verbose=True))
