##############################################################
# Julia type for Mori dream spaces that are toric varieties
# #############################################################
abstract type ToricVarietyMDS <: MoriDreamSpace end

cox_ring_relations(::ToricVarietyMDS) = []


