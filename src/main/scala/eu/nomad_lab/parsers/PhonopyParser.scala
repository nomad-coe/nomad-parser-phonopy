/*
   Copyright 2016-2017 The NOMAD Developers Group

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
 */
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
          (lab.PhonopyVersionInfo.toMap.map {
            case (key, value) =>
              (key -> jn.JString(value.toString))
          }(breakOut): List[(String, jn.JString)])
      )) :: Nil
  ),
  mainFileTypes = Seq("text/plain"),
  ancillaryFilesPrefilter = AncillaryFilesPrefilter.ParentSubtree,
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
  ) ++ lab.DefaultPythonInterpreter.commonDirMapping(),
  metaInfoEnv = Some(lab.meta.KnownMetaInfoEnvs.phonopy)
) {
  val fileRe = ".*/phonopy-FHI-aims-displacement-0*1/control\\.in$".r
  override def isMainFile(filePath: String, bytePrefix: Array[Byte], stringPrefix: Option[String]): Option[ParserMatch] = {
    filePath match {
      case fileRe() =>
        Some(ParserMatch(mainFileMatchPriority, mainFileMatchWeak))
      case _ =>
        None
    }
  }
}
