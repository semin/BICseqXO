import sbt._

object PluginDef extends Build {
    override lazy val projects = Seq(root)
    lazy val root = Project("plugins", file(".")) dependsOn(webPlugin)
    lazy val webPlugin = uri("git://github.com/eed3si9n/sbt-assembly.git")
}
// vim: set ts=4 sw=4 et:
