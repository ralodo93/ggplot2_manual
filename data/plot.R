top_10_up <-
  results %>% arrange(PValue) %>% filter(logFC > 1) %>% head(5)
top_10_down <-
  results %>% arrange(PValue) %>% filter(logFC < 1) %>% head(5)

top_10_up_top_10_down <- rbind(top_10_up, top_10_down)

expression_values <-
  expression_matrix %>% filter(gene %in% top_10_up_top_10_down$gene) %>%
  left_join(sampleinfo) %>% 
  mutate(gene = factor(gene, levels = top_10_up_top_10_down$gene)) %>%
  mutate(SampleName = factor(sampleinfo %>% arrange(CellType) %>% pull(SampleName)))


ggplot(expression_values, aes(y = gene, x = SampleName, fill = expression)) +
  theme_minimal()+
  geom_tile(color = "gray50")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
