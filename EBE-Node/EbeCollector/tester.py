
import EbeCollector

reader = EbeCollector.EbeDBReader("testDB/collected.db")

res = reader.getAttendance()
print(res)
