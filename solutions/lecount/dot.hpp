#ifndef DOT_H
#define DOT_H

struct DOTAttribute
{
	std::string name;
	std::string value;
	
	DOTAttribute(const std::string &name, const std::string &value) : name(name), value(value) {}
	
	void print_to_file(FILE *f)
	{
		// NOTE: not escaped
		fprintf(f, "%s=\"%s\"", name.c_str(), value.c_str());
	}
};

struct DOTElement
{
	std::vector<DOTAttribute*> attributes;
	
	~DOTElement()
	{
		for (unsigned i = 0; i < attributes.size(); ++i) {
			delete attributes[i];
		}
	}
	
	void add_attribute(const std::string &name, const std::string &value)
	{
		DOTAttribute *attribute = new DOTAttribute(name, value);
		attributes.push_back(attribute);
	}
	
	void print_attributes_to_file(FILE *f)
	{
		if (attributes.size() > 0) {
			fprintf(f, " [");
			for (unsigned i = 0; i < attributes.size(); ++i) {
				if (i > 0) fprintf(f, ",");
				attributes[i]->print_to_file(f);
			}
			fprintf(f, "]");
		}
	}
};

struct DOTNode : public DOTElement
{
	std::string name;
	
	DOTNode(const std::string &name) : name(name) {}
	
	void print_to_file(FILE *f)
	{
		fprintf(f, "%s", name.c_str());
		print_attributes_to_file(f);
		fprintf(f, ";\n");
	}
};

struct DOTEdge : public DOTElement
{
	std::string start;
	std::string end;
	
	DOTEdge(const std::string &start, const std::string &end) : start(start), end(end) {}
	
	void print_to_file(FILE *f)
	{
		fprintf(f, "%s -> %s", start.c_str(), end.c_str());
		print_attributes_to_file(f);
		fprintf(f, ";\n");
	}
};

struct DOTDigraph
{
	std::vector<DOTNode*> nodes;
	std::vector<DOTEdge*> edges;
	
	DOTDigraph() {}
	
	~DOTDigraph()
	{
		for (unsigned i = 0; i < nodes.size(); ++i) {
			delete nodes[i];
		}
		
		for (unsigned i = 0; i < edges.size(); ++i) {
			delete edges[i];
		}
	}
	
	DOTNode *add_node(const std::string &name)
	{
		DOTNode *node = new DOTNode(name);
		nodes.push_back(node);
		return node;
	}
	
	DOTEdge *add_edge(const std::string &start, const std::string &end)
	{
		DOTEdge *edge = new DOTEdge(start, end);
		edges.push_back(edge);
		return edge;
	}
	
	void print_to_file(FILE *f)
	{
		fprintf(f, "digraph G {\n");
		
		for (unsigned i = 0; i < nodes.size(); ++i) {
			nodes[i]->print_to_file(f);
		}
		
		for (unsigned i = 0; i < edges.size(); ++i) {
			edges[i]->print_to_file(f);
		}
		
		fprintf(f, "}\n");
	}
};

#endif
