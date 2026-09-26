// File: Transfer/src/DynamoEngine/UI/Widgets/UILabel.cpp

#include "DynamoEngine/UI/Widgets/UILabel.hpp"

namespace DynamoEngine
{
void UILabel::update(float delta_seconds)
{
    if (m_text_source) // a text source is set: refresh the text every frame
    {
        m_text = m_text_source();
    }
}

void UILabel::draw(UIGeometryBuilder& builder) const
{
    builder.addText(m_text, {m_rect.x, m_rect.y});
}
} // namespace DynamoEngine
