
import EbeCollector

reader = EbeCollector.EbeDBReader("testDB/collected.db")

res = reader.evaluateExpression("Psi_2(pion)")
print(res)
