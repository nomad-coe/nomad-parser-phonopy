package eu.nomad_lab.parsers

import org.specs2.mutable.Specification

object PhonopyParserSpec extends Specification {
  "PhonopyParserTest" >> {
    "Si test with json-events" >> {
      ParserRun.parse(PhonopyParser, "parsers/phonopy/test/examples/phonopy-FHI-aims-displacement-01/control.in", "json-events") must_== ParseResult.ParseSuccess
    }
    "Si test with json" >> {
      ParserRun.parse(PhonopyParser, "parsers/phonopy/test/examples/phonopy-FHI-aims-displacement-01/control.in", "json") must_== ParseResult.ParseSuccess
    }
  }
}
