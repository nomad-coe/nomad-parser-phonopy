package eu.nomad_lab.parsers
import eu.{ nomad_lab => lab }
import org.{ json4s => jn }
import scala.collection.breakOut

object PhonopyParser extends SimpleExternalParserGenerator(
  name = "PhonopyParser",
  parserInfo = jn.JObject(
    ("name" -> jn.JString("PhonopyParser")) ::
      ("parserId" -> jn.JString("PhonopyParser" + lab.PhonopyVersionInfo.version)) ::
      ("versionInfo" -> jn.JObject(
        ("nomadCoreVersion" -> jn.JObject(lab.NomadCoreVersionInfo.toMap.map {
          case (k, v) => k -> jn.JString(v.toString)
        }(breakOut): List[(String, jn.JString)])) ::
          (lab.PhonophyVersionInfo.toMap.map {
            case (key, value) =>
              (key -> jn.JString(value.toString))
          }(breakOut): List[(String, jn.JString)])
      )) :: Nil
  ),
  mainFileTypes = Seq("application/x-gtar"),
  mainFileRe = "".r,
  cmd = Seq(lab.DefaultPythonInterpreter.pythonExe(), "${envDir}/parsers/phonopy/parser/parser-phonopy/Get_Force_Constants.py", "${mainFileUri}", "${mainFilePath}"),
  resList = Seq(
    "parser-phonopy/con.py",
    "parser-phonopy/Get_Force_Constants.py",
    "parser-phonopy/Get_Properties.py",
    "parser-phonopy/PhononModulesNomad.py",
    "parser-phonopy/setup_paths.py",
    "nomad_meta_info/public.nomadmetainfo.json",
    "nomad_meta_info/common.nomadmetainfo.json",
    "nomad_meta_info/meta_types.nomadmetainfo.json",
    "nomad_meta_info/phonopy.nomadmetainfo.json"
  ) ++ lab.DefaultPythonInterpreter.commonFiles(),
  dirMap = Map(
    "parser-phonopy" -> "parsers/phonopy/parser/parser-phonopy",
    "nomad_meta_info" -> "nomad-meta-info/meta_info/nomad_meta_info"
  ) ++ lab.DefaultPythonInterpreter.commonDirMapping()
) {
  override def isMainFile(filePath: String, bytePrefix: Array[Byte], stringPrefix: Option[String]): Option[ParserMatch] = {
    if (filePath.endsWith(".*/phonopy-FHI-aims-displacement-0*1/control\\.in$"))
      Some(ParserMatch(mainFileMatchPriority, mainFileMatchWeak))
    else
      None
  }
}
